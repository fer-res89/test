"""The `process`, `pipe-all`, and `run-dirs` subcommands: run MCMCTree +
CodeML on one FASTA, on every FASTA in the current directory, or across
subdirectories.

Every run happens in its own scratch directory (``_aeon_work_<stem>``), so
MCMCTree/CodeML's fixed output names (tmp0001.*, rst2, in.BV, mcmc.txt,
out.txt, FigTree.tre, wag.dat) can never collide between genes, and runs can
go in parallel (``--jobs``) even when many FASTAs share a directory. Only the
renamed outputs are moved to the destination:

    <stem>.FigTree.tre  <stem>.mcmc.txt  <stem>.out.txt  <stem>.ctl  [<stem>.in.BV]

Prior-only runs (``--prior-only``: usedata = 0, no CodeML step) are named
``<stem>.prior.*`` unless ``--out-dir`` sends them elsewhere.
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from ._common import (_load_ctl_template, _write_wag, ctl_get, ctl_set, format_root_age,
                      generate_ctl, run_cmd)
from .sanitize_tree import parse_sanitize_letters, sanitize_tree_file

OUTPUTS = (("FigTree.tre", ".FigTree.tre"), ("mcmc.txt", ".mcmc.txt"), ("out.txt", ".out.txt"))


# ---------------------------------------------------------------------------
# Tree helpers
# ---------------------------------------------------------------------------

def _first_tree(text: str) -> str:
    """First Newick tree in a file, dropping an MCMCTree '<ntax> <ntree>' header."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if lines and re.fullmatch(r"\s*\d+\s+\d+\s*", lines[0]):
        lines = lines[1:]
    body = " ".join(lines)
    q = False
    depth = 0
    for i, ch in enumerate(body):
        if q:
            q = ch != "'"
        elif ch == "'":
            q = True
        elif ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
        elif ch == ";" and depth == 0:
            return body[: i + 1]
    return body.rstrip() + ";"


def calibrated_tree_text(tree_text: str, root_cal: str) -> str:
    """MCMCTree tree file: header, topology without branch lengths or unquoted
    (support) labels, quoted calibrations kept, `root_cal` placed on the root."""
    t = _first_tree(tree_text)
    t = re.sub(r"\[[^\]]*\]", "", t)                       # comments
    t = re.sub(r":\s*[-+0-9.eE]+", "", t)                  # branch lengths
    t = re.sub(r"\)\s*(?!')[^,();']+", ")", t)             # unquoted internal labels
    t = t.rstrip().rstrip(";").rstrip()
    # replace whatever follows the root's closing parenthesis
    t = re.sub(r"\)\s*('[^']*')?\s*$", ")", t)
    cal = root_cal.strip().strip("'")
    ntax = re.sub(r"'[^']*'", "", t).count(",") + 1        # commas outside quoted labels
    return f"{ntax} 1\n{t}'{cal}';\n"


def ensure_header(tree_text: str) -> str:
    """MCMCTree needs a '<ntaxa> <ntrees>' first line; add one if missing."""
    lines = [ln for ln in tree_text.splitlines() if ln.strip()]
    if lines and re.fullmatch(r"\s*\d+\s+\d+\s*", lines[0]):
        return tree_text
    t = _first_tree(tree_text)
    ntax = re.sub(r"'[^']*'|\[[^\]]*\]", "", t).count(",") + 1
    return f"{ntax} 1\n{t}\n"


def _root_is_calibrated(tree_text: str) -> bool:
    t = _first_tree(tree_text).rstrip(";").rstrip()
    return bool(re.search(r"\)\s*('[^']*'|[<>][^,();]*)\s*$", t))


# ---------------------------------------------------------------------------
# One gene
# ---------------------------------------------------------------------------

def do_process(seqfile: str, treefile: str, ctl_template: str | None, wag_file: str | None,
               *, prior_only: bool = False, root_age: str | None = None, root_cal: str | None = None,
               out_dir: str | None = None, overrides: dict | None = None, keep_work: bool = False) -> bool:
    """Run one gene. Returns True on success. Failed runs keep their work dir."""
    seq = Path(seqfile).resolve()
    tree = Path(treefile).resolve()
    if not seq.exists() or not tree.exists():
        print(f"[ERROR] sequence or tree not found: {seq} / {tree}")
        return False
    dest = Path(out_dir).resolve() if out_dir else seq.parent
    dest.mkdir(parents=True, exist_ok=True)
    stem = seq.stem + (".prior" if prior_only and not out_dir else "")
    work = dest / f"_aeon_work_{stem}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir()

    # inputs
    shutil.copy(seq, work / seq.name)
    tree_text = tree.read_text()
    tree_text = calibrated_tree_text(tree_text, root_cal) if root_cal else ensure_header(tree_text)
    (work / "tree.tre").write_text(tree_text)

    template_text = _load_ctl_template(ctl_template)
    ov = {"mcmcfile": "mcmc.txt", "outfile": "out.txt",
          "usedata": "0    * prior only" if prior_only else "3"}
    if root_age is not None:
        ov["RootAge"] = None if root_age.lower() == "none" else format_root_age(root_age)
    elif root_cal and ctl_get(template_text, "RootAge"):
        ov["RootAge"] = None
        print(f"[INFO] {seq.name}: --root-cal given, RootAge removed from the control file")
    ov.update(overrides or {})
    ctl = work / "run.ctl"
    ctl_text = generate_ctl(template_text, seq.name, "tree.tre", ctl, ov)
    if not ctl_get(ctl_text, "RootAge") and not _root_is_calibrated(tree_text):
        print(f"[WARN] {seq.name}: no RootAge and no calibration on the root; "
              f"MCMCTree will stop (use --root-cal or --root-age)")
    _write_wag(work, wag_file)

    def fail(msg):
        print(f"[ERROR] {seq.name}: {msg}; work dir kept: {work}")
        return False

    if prior_only:
        try:
            run_cmd(["mcmctree", "run.ctl"], cwd=str(work))
        except subprocess.CalledProcessError:
            return fail("mcmctree (prior only) failed")
    else:
        # 1) usedata = 3: MCMCTree writes tmp0001.* for CodeML
        try:
            run_cmd(["mcmctree", "run.ctl"], cwd=str(work))
        except subprocess.CalledProcessError:
            files = [work / f for f in ("tmp0001.ctl", "tmp0001.txt", "tmp0001.trees")]
            if not all(f.exists() and f.stat().st_size > 0 for f in files):
                return fail("first mcmctree failed and tmp0001.* missing")
            print(f"[WARNING] first mcmctree failed on {seq.name}, outputs present -> continuing.")
        tmp = work / "tmp0001.ctl"
        if not tmp.exists():
            return fail("tmp0001.ctl missing")
        txt = tmp.read_text()
        txt = re.sub(r"(?m)^\s*(aaRatefile|fix_alpha|alpha|ncatG)\s*=.*\n?", "", txt)
        txt = re.sub(r"(?m)^\s*model\s*=.*$",
                     "model = 2\naaRatefile = wag.dat\nfix_alpha = 0\nalpha = .5\nncatG = 4", txt)
        tmp.write_text(txt)
        # 2) CodeML: gradient and Hessian -> rst2
        rst2 = work / "rst2"
        try:
            run_cmd(["codeml", "tmp0001.ctl"], cwd=str(work))
        except subprocess.CalledProcessError:
            if not (rst2.exists() and rst2.stat().st_size > 0):
                return fail("codeml failed and rst2 missing")
            print(f"[WARNING] codeml failed on {seq.name}, but rst2 present -> continuing.")
        if not rst2.exists():
            return fail("rst2 missing after codeml")
        rst2.rename(work / "in.BV")
        # 3) usedata = 2: MCMC with the approximate likelihood
        ctl.write_text(ctl_set(ctl.read_text(), "usedata", "2 in.BV"))
        try:
            run_cmd(["mcmctree", "run.ctl"], cwd=str(work))
        except subprocess.CalledProcessError:
            return fail("final mcmctree failed")

    if not (work / "FigTree.tre").exists():
        return fail("FigTree.tre not produced")
    for src, ext in OUTPUTS:
        p = work / src
        if p.exists():
            shutil.move(str(p), dest / f"{stem}{ext}")
    shutil.copy(ctl, dest / f"{stem}.ctl")
    if (work / "in.BV").exists():
        shutil.move(str(work / "in.BV"), dest / f"{stem}.in.BV")
    if not keep_work:
        shutil.rmtree(work, ignore_errors=True)
    return True


def _expected_fig(fa: Path, prior_only: bool, out_dir: str | None) -> Path:
    dest = Path(out_dir) if out_dir else fa.parent
    return dest / f"{fa.stem}{'.prior' if prior_only and not out_dir else ''}.FigTree.tre"


def _run_opts(args) -> dict:
    ov = {}
    for key, attr in (("burnin", "burnin"), ("sampfreq", "sampfreq"), ("nsample", "nsample"),
                      ("seed", "seed")):
        v = getattr(args, attr, None)
        if v is not None:
            ov[key] = str(v)
    return dict(prior_only=args.prior_only, root_age=args.root_age, root_cal=args.root_cal,
                overrides=ov, keep_work=args.keep_work)


def _worker(task):
    fa, tree, ctl_t, wag, out_dir, opts = task
    try:
        return fa, do_process(fa, tree, ctl_t, wag, out_dir=out_dir, **opts)
    except Exception as e:  # never kill the pool
        print(f"[ERROR] {fa}: {e}")
        return fa, False


def _run_tasks(tasks, jobs: int):
    ok = bad = 0
    failed = []
    if jobs <= 1:
        results = map(_worker, tasks)
        for fa, res in results:
            ok += res; bad += (not res)
            if not res:
                failed.append(fa)
    else:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for fut in as_completed([ex.submit(_worker, t) for t in tasks]):
                fa, res = fut.result()
                ok += res; bad += (not res)
                if not res:
                    failed.append(fa)
    print(f"[SUMMARY] {ok} succeeded, {bad} failed")
    for f in failed:
        print("  failed:", f)
    return bad == 0


def _resolve_tree_for_fasta(fa_path: Path, explicit_treefile: str | None, tree_dir: str | None,
                            tree_ext: str = '.tre') -> Path | None:
    """Explicit tree file first, else <fa.stem><tree_ext> in tree_dir."""
    if explicit_treefile:
        p = Path(explicit_treefile)
        if p.is_file():
            return p
        print(f"[WARN] Explicit tree file not found: {p}")
    if tree_dir:
        cand = Path(tree_dir) / f"{fa_path.stem}{tree_ext}"
        if cand.is_file():
            return cand
        print(f"[WARN] No matching tree in --tree-dir for {fa_path.name}: expected {cand.name}")
    return None


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_process(args):
    fa = Path(args.seqfile)
    tree_path = _resolve_tree_for_fasta(fa, args.treefile, args.tree_dir, args.tree_ext)
    if not tree_path:
        sys.exit("[ERROR] Could not resolve a tree file (provide a treefile or --tree-dir).")
    return do_process(str(fa), str(tree_path), args.ctl_template, args.wag_file,
                      out_dir=args.out_dir, **_run_opts(args))


def cmd_pipe_all(args):
    cwd = Path.cwd()
    tasks = []
    for fa in sorted(cwd.glob(f'*{args.fa_ext}')):
        if args.resume and _expected_fig(fa, args.prior_only, args.out_dir).exists():
            continue
        tree_path = _resolve_tree_for_fasta(fa, args.treefile, args.tree_dir, args.tree_ext)
        if not tree_path:
            print(f"[WARN] Skipping {fa.name}: no tree found (need a treefile or --tree-dir match).")
            continue
        tasks.append((str(fa), str(tree_path), args.ctl_template, args.wag_file, args.out_dir,
                      _run_opts(args)))
    return _run_tasks(tasks, args.jobs)


def cmd_run_dirs(args):
    cwd = Path.cwd()
    missing = []
    tasks = []
    opts = _run_opts(args)
    for sub in sorted(cwd.iterdir()):
        if not sub.is_dir() or sub.name.startswith("_aeon_work_"):
            continue
        out_dir = str(Path(args.out_dir) / sub.name) if args.out_dir else None
        if args.tree_dir:
            tree_dir = Path(args.tree_dir)
            if not tree_dir.is_dir():
                sys.exit(f"[ERROR] --tree-dir not found: {tree_dir}")
            search_root = tree_dir / sub.name if args.tree_dir_search else tree_dir
            if args.tree_dir_search and not search_root.is_dir():
                msg = f"{sub.name}: expected tree dir {search_root} (from --tree-dir-search), but not found"
                print(f"[ERROR] {msg}"); missing.append(msg); continue
            for fa in sorted(sub.glob(f"*{args.fa_ext}")):
                if args.resume and _expected_fig(fa, args.prior_only, out_dir).exists():
                    continue
                cand = search_root / f"{fa.stem}{args.tree_ext}"
                if not cand.is_file():
                    msg = f"{sub.name}/{fa.name}: expected {cand.name} in {search_root}, but not found"
                    print(f"[ERROR] {msg}"); missing.append(msg); continue
                tree_to_use = cand
                if args.sanitize is not None:
                    san = parse_sanitize_letters(args.sanitize)
                    sanitized = fa.parent / (cand.stem + ".sanitized.tre")
                    sanitize_tree_file(cand, sanitized, collapse=san["collapse"],
                                       keep_support=san["keep_support"], keep_lengths=san["keep_lengths"],
                                       allow_polytomies=san["allow_polytomies"], add_header=san["add_header"],
                                       quiet=san["quiet"])
                    tree_to_use = sanitized
                tasks.append((str(fa), str(tree_to_use), args.ctl_template, args.wag_file, out_dir, opts))
        else:
            top_tree = cwd / f"{sub.name}.tre"
            if top_tree.exists():
                shutil.copy(top_tree, sub / top_tree.name)
            for fa in sorted(sub.glob(f"*{args.fa_ext}")):
                if args.resume and _expected_fig(fa, args.prior_only, out_dir).exists():
                    continue
                if top_tree.exists():
                    tasks.append((str(fa), str(sub / top_tree.name), args.ctl_template, args.wag_file,
                                  out_dir, opts))
                else:
                    print(f"[WARN] {sub.name}/{fa.name}: no <subdir>.tre and no --tree-dir; skipping.")
                    missing.append(f"{sub.name}/{fa.name}: no <subdir>.tre")
    ok = _run_tasks(tasks, args.jobs)
    if missing:
        print("\n[SUMMARY] Missing tree matches:")
        for m in missing:
            print(" -", m)
    return ok


# ---------------------------------------------------------------------------
# Argument registration
# ---------------------------------------------------------------------------

def _add_run_options(p, parallel: bool):
    g = p.add_argument_group("run options")
    g.add_argument('--ctl-template', help='Use this .ctl template instead of the embedded default')
    g.add_argument('--wag-file', help='Amino-acid matrix to use instead of the embedded WAG '
                                      '(any PAML matrix; it is copied as wag.dat)')
    g.add_argument('--prior-only', action='store_true',
                   help='Sample from the prior (usedata = 0, no CodeML). Outputs are named '
                        '<stem>.prior.* unless --out-dir is given')
    g.add_argument('--root-age', metavar='VALUE',
                   help="Set RootAge in the control file, e.g. '<3.0' or '>2.6<2.9'; "
                        "'none' removes it. Default: as in the template (embedded: absent)")
    g.add_argument('--root-cal', metavar='CAL',
                   help="Put a calibration on the root of the tree, e.g. 'B(2.6,2.9,0.025,0.025)' "
                        "(MCMCTree time units). Writes header, strips lengths and support "
                        "values, keeps other quoted calibrations, and removes RootAge "
                        "from the control file unless --root-age is also given")
    g.add_argument('--out-dir', help='Write outputs here instead of next to the FASTA '
                                     '(run-dirs: <out-dir>/<subdir>/)')
    g.add_argument('--burnin', type=int, help='Override burnin')
    g.add_argument('--sampfreq', type=int, help='Override sampfreq')
    g.add_argument('--nsample', type=int, help='Override nsample')
    g.add_argument('--seed', type=int, help='Override seed (e.g. for a second chain); template default -1')
    g.add_argument('--keep-work', action='store_true', help='Keep the _aeon_work_<stem> scratch directory')
    if parallel:
        g.add_argument('-j', '--jobs', type=int, default=1, help='Genes to run in parallel (default 1)')


def register(sub) -> None:
    """Register the `process`, `pipe-all`, and `run-dirs` subcommands."""
    p1 = sub.add_parser('process', help='Run MCMCTree + CodeML on a single file')
    p1.add_argument('seqfile', help='Input FASTA file')
    p1.add_argument('treefile', nargs='?', help='Tree file (optional if --tree-dir is provided)')
    p1.add_argument('--tree-dir', help='Directory containing per-gene trees; match by FASTA basename')
    p1.add_argument('--tree-ext', default='.tre', help='Tree filename extension in --tree-dir (default: .tre)')
    _add_run_options(p1, parallel=False)
    p1.set_defaults(func=cmd_process)

    p2 = sub.add_parser('pipe-all', help='Run pipeline on all FASTAs in cwd')
    p2.add_argument('treefile', nargs='?', help='Tree file for all FASTAs (optional if --tree-dir is used)')
    p2.add_argument('-e', '--fa-ext', default='.fa', help='FASTA extension (default: .fa)')
    p2.add_argument('-r', '--resume', action='store_true', help='Skip FASTAs whose FigTree.tre output exists')
    p2.add_argument('--tree-dir', help='Directory containing per-gene trees; match by FASTA basename')
    p2.add_argument('--tree-ext', default='.tre', help='Tree filename extension in --tree-dir (default: .tre)')
    _add_run_options(p2, parallel=True)
    p2.set_defaults(func=cmd_pipe_all)

    p3 = sub.add_parser(
        'run-dirs',
        help='Run pipeline across subdirectories',
        description=(
            "Run the pipeline across all subdirectories in the current working directory.\n\n"
            "For each subdirectory:\n"
            "  • The program looks for a tree file in the parent directory with the same base name "
            "as the subdirectory (e.g. subdir 'pal_neog' expects 'pal_neog.tre').\n"
            "  • If such a file exists, it is copied into the subdirectory and used for all FASTAs there.\n\n"
            "  • Use --tree-dir (and --tree-ext if necessary) to use a tree for each FASTA by matching "
            "basename from a directory.\n"
            "       example: for FASTA OG0001.fa and --tree-dir trees/ --tree-ext .tre, we look for trees/OG0001.tre.\n"
            "  • With --tree-dir-search, trees for subdirectory X are searched under PATH(--tree-dir)/X/.\n\n"
            "Each gene runs in its own scratch directory, so --jobs N is safe."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p3.add_argument('-e', '--fa-ext', default='.fa', help='FASTA extension (default: .fa)')
    p3.add_argument('-r', '--resume', action='store_true', help='Skip FASTAs whose FigTree.tre output exists')
    p3.add_argument('--tree-dir', help='Directory containing per-gene trees; match by FASTA basename')
    p3.add_argument('--tree-ext', default='.tre', help='Tree filename extension in --tree-dir (default: .tre)')
    p3.add_argument('--sanitize', nargs='?', const='',
                    help=("Sanitize per-gene trees from --tree-dir using letters "
                          "(l=keep-lengths, s=keep-support, p=allow-polytomies, h=add-header, c=no-collapse). "
                          "Example: '--sanitize sc'"))
    p3.add_argument("--tree-dir-search", action="store_true",
                    help=("Search for trees under PATH(--tree-dir)/<subdir_name>/ when matching trees by "
                          "FASTA basename. Requires --tree-dir."))
    _add_run_options(p3, parallel=True)
    p3.set_defaults(func=cmd_run_dirs)
