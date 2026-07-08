"""The `process`, `pipe-all`, and `run-dirs` subcommands: run MCMCTree +
CodeML on one FASTA, on every FASTA in the current directory, or across
subdirectories.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

from ._common import run_cmd, generate_ctl, _load_ctl_template, _write_wag
from .sanitize_tree import parse_sanitize_letters, sanitize_tree_file

def do_process(seqfile: str, treefile: str, ctl_template: str | None, wag_file: str | None):
    seq = Path(seqfile)
    tree = Path(treefile)
    if not seq.exists() or not tree.exists():
        sys.exit("ERROR: sequence or tree not found")
    dest = seq.parent
    dest_ctl = dest / seq.with_suffix('.ctl').name

    template_text = _load_ctl_template(ctl_template)
    generate_ctl(template_text, seq.name, str(tree.resolve()), dest_ctl)


    _write_wag(dest, wag_file)

    # First MCMCTree run (allow failure if outputs are present)
    try:
        run_cmd(['mcmctree', str(dest_ctl)], cwd=str(dest))
    except subprocess.CalledProcessError:
        files = [dest / f for f in ("tmp0001.ctl", "tmp0001.txt", "tmp0001.trees", "tmp0001.out")]
        if all(f.exists() and f.stat().st_size > 0 for f in files):
            print(f"[WARNING] first mcmctree failed on {seq.name}, outputs present → continuing.")
        else:
            print(f"[ERROR] first mcmctree failed and outputs missing → skipping.")
            return

    # Rename tmp0001.ctl to <seq>.tmp0001.ctl and patch for CodeML
    tmp = dest / 'tmp0001.ctl'
    newtmp = dest / seq.with_suffix('.tmp0001.ctl').name
    if tmp.exists():
        tmp.rename(newtmp)
    txt = newtmp.read_text()
    # Ensure we point to wag.dat that we wrote/copied above
    txt = re.sub(r'(?m)^\s*model\s*=.*$',
                 'model = 2\naaRatefile = wag.dat\nfix_alpha = 0\nalpha = .5\nncatG = 4', txt)
    newtmp.write_text(txt)

    # Run CodeML (allow continuation if rst2 exists)
    rst2 = dest / 'rst2'
    try:
        run_cmd(['codeml', str(newtmp)], cwd=str(dest))
    except subprocess.CalledProcessError:
        if rst2.exists() and rst2.stat().st_size > 0:
            print(f"[WARNING] codeml failed on {seq.name}, but rst2 present → continuing.")
        else:
            print(f"[ERROR] codeml failed and rst2 missing → skipping.")
            return

    if rst2.exists():
        rst2.rename(dest / 'in.BV')

    # Flip usedata=3 → 2 to consume in.BV
    ctl_text = dest_ctl.read_text()
    ctl_text = re.sub(r'(?m)^\s*usedata\s*=\s*3', 'usedata = 2', ctl_text)
    dest_ctl.write_text(ctl_text)

    # Final MCMCTree
    run_cmd(['mcmctree', str(dest_ctl)], cwd=str(dest))

    # Rename outputs
    for src, ext in [('mcmc.txt', '.mcmc.txt'), ('out.txt', '.out.txt'), ('FigTree.tre', '.FigTree.tre')]:
        p = dest / src
        if p.exists():
            p.rename(dest / seq.with_suffix(ext).name)
    ib = dest / 'in.BV'
    if ib.exists():
        ib.rename(dest / seq.with_suffix('.in.BV').name)

    # Cleanup
    for fname in ["tmp0001.txt", "tmp0001.trees", "tmp0001.out",
                  "SeedUsed", "lnf", "out.BV", "rates", "rst", "rst1", "rub"]:
        f = dest / fname
        if f.exists():
            try:
                f.unlink()
            except Exception as e:
                print(f"[WARNING] could not delete {fname}: {e}")

def _resolve_tree_for_fasta(fa_path: Path,
                            explicit_treefile: str | None,
                            tree_dir: str | None,
                            tree_ext: str = '.tre') -> Path | None:
    """
    Decide which tree file to use for a given FASTA:
      Precedence:
        1) explicit_treefile (if provided)
        2) match by basename in tree_dir:  <fa.stem><tree_ext>
      Returns a Path or None if not found.
    """
    # 1) explicit
    if explicit_treefile:
        p = Path(explicit_treefile)
        if p.is_file():
            return p
        print(f"[WARN] Explicit tree file not found: {p}")

    # 2) match by basename from --tree-dir
    if tree_dir:
        cand = Path(tree_dir) / f"{fa_path.stem}{tree_ext}"
        if cand.is_file():
            return cand
        print(f"[WARN] No matching tree in --tree-dir for {fa_path.name}: expected {cand.name}")

    return None

def cmd_process(args):
    fa = Path(args.seqfile)
    tree_path = _resolve_tree_for_fasta(fa, args.treefile, args.tree_dir, args.tree_ext)
    if not tree_path:
        sys.exit("[ERROR] Could not resolve a tree file (provide --treefile or --tree-dir).")
    do_process(str(fa), str(tree_path), args.ctl_template, args.wag_file)



def cmd_pipe_all(args):
    cwd = Path.cwd()
    for fa in cwd.glob(f'*{args.fa_ext}'):
        fig = cwd / f"{fa.stem}.FigTree.tre"
        if args.resume and fig.exists():
            continue
        tree_path = _resolve_tree_for_fasta(fa, args.treefile, args.tree_dir, args.tree_ext)
        if not tree_path:
            print(f"[WARN] Skipping {fa.name}: no tree found (need --treefile or --tree-dir match).")
            continue
        do_process(str(fa), str(tree_path), args.ctl_template, args.wag_file)



def cmd_run_dirs(args):
    cwd = Path.cwd()
    missing = []  # collect missing matches for a summary

    if args.tree_dir:
        # STRICT per-gene mode (no <subdir>.tre fallback)
        tree_dir = Path(args.tree_dir)
        if not tree_dir.is_dir():
            sys.exit(f"[ERROR] --tree-dir not found: {tree_dir}")

        for sub in cwd.iterdir():
            if not sub.is_dir():
                continue

            # Determine where to search for trees for THIS subdir
            # Default: the root --tree-dir
            search_root = tree_dir

            # If requested, scope search to --tree-dir/<subdir_name>/
            if args.tree_dir_search:
                scoped = tree_dir / sub.name
                if not scoped.is_dir():
                    msg = (f"{sub.name}: expected tree dir {scoped} (from --tree-dir-search), but not found")
                    print(f"[ERROR] {msg}")
                    missing.append(msg)
                    continue
                search_root = scoped

            for fa in sub.glob(f"*{args.fa_ext}"):
                fig = sub / f"{fa.stem}.FigTree.tre"
                if args.resume and fig.exists():
                    continue

                cand = search_root / f"{fa.stem}{args.tree_ext}"
                if not cand.is_file():
                    msg = (f"{sub.name}/{fa.name}: expected {cand.name} in {search_root}, but not found")
                    print(f"[ERROR] {msg}")
                    missing.append(msg)
                    continue

                tree_to_use = cand
                if args.sanitize is not None:
                    san = parse_sanitize_letters(args.sanitize)
                    # Write sanitized copy next to the FASTA (do not mutate originals)
                    sanitized = fa.parent / (cand.stem + ".sanitized.tre")
                    sanitize_tree_file(
                        cand, sanitized,
                        collapse=san["collapse"],
                        keep_support=san["keep_support"],
                        keep_lengths=san["keep_lengths"],
                        allow_polytomies=san["allow_polytomies"],
                        add_header=san["add_header"],
                        quiet=san["quiet"]
                    )
                    tree_to_use = sanitized

                do_process(str(fa), str(tree_to_use), args.ctl_template, args.wag_file)

    else:
        # ORIGINAL mode: use <subdir>.tre from parent if present
        for sub in cwd.iterdir():
            if not sub.is_dir():
                continue

            top_tree = cwd / f"{sub.name}.tre"
            have_top = top_tree.exists()
            if have_top:
                shutil.copy(top_tree, sub / top_tree.name)

            for fa in sub.glob(f"*{args.fa_ext}"):
                fig = sub / f"{fa.stem}.FigTree.tre"
                if args.resume and fig.exists():
                    continue

                if have_top:
                    tree_path = sub / top_tree.name
                    do_process(str(fa), str(tree_path), args.ctl_template, args.wag_file)
                else:
                    print(f"[WARN] {sub.name}/{fa.name}: no <subdir>.tre and no --tree-dir; skipping.")
                    missing.append(f"{sub.name}/{fa.name}: no <subdir>.tre")

    # Summary for visibility
    if missing:
        print("\n[SUMMARY] Missing tree matches:")
        for m in missing:
            print(" -", m)


# =========================
# HPD extraction
# =========================



def register(sub) -> None:
    """Register the `process`, `pipe-all`, and `run-dirs` subcommands."""
    # process
    p1 = sub.add_parser('process', help='Run MCMCTree + CodeML on a single file')
    p1.add_argument('seqfile', help='Input FASTA file')
    p1.add_argument('treefile', nargs='?', help='Tree file (optional if --tree-dir is provided)')
    p1.add_argument('--ctl-template', help='Use this .ctl template instead of embedded default')
    p1.add_argument('--wag-file', help='Use this WAG matrix file instead of embedded default')
    p1.add_argument('--tree-dir', help='Directory containing per-gene trees; match by FASTA basename')
    p1.add_argument('--tree-ext', default='.tre', help='Tree filename extension in --tree-dir (default: .tre)')
    p1.set_defaults(func=cmd_process)

    # pipe-all
    p2 = sub.add_parser('pipe-all', help='Run pipeline on all FASTAs in cwd')
    p2.add_argument('treefile', nargs='?', help='Tree file for all FASTAs (optional if --tree-dir is used)')
    p2.add_argument('-e', '--fa-ext', default='.fa', help='FASTA extension (default: .fa)')
    p2.add_argument('-r', '--resume', action='store_true', help='Skip FASTAs with existing .FigTree.tre')
    p2.add_argument('--ctl-template', help='Override .ctl template')
    p2.add_argument('--wag-file', help='Override WAG matrix file')
    p2.add_argument('--tree-dir', help='Directory containing per-gene trees; match by FASTA basename')
    p2.add_argument('--tree-ext', default='.tre', help='Tree filename extension in --tree-dir (default: .tre)')
    p2.set_defaults(func=cmd_pipe_all)

    # run-dirs
    p3 = sub.add_parser(
    'run-dirs',
    help='Run pipeline across subdirectories',
    description=(
        "Run the pipeline across all subdirectories in the current working directory.\n\n"
        "For each subdirectory:\n"
        "  • The program looks for a tree file in the parent directory with the same base name "
        "as the subdirectory (e.g. subdir 'pal_neog' expects 'pal_neog.tre').\n"
        "  • If such a file exists, it is copied into the subdirectory and used for all FASTAs there.\n"
        "  • Only subdirectories with a matching <subdir>.tre will have a tree copied in.\n\n"
        "Example:\n"
        "  parent_dir/\n"
        "    ├── pal_neog/\n"
        "    │     ├── OG0001.fa\n"
        "    │     └── OG0002.fa\n"
        "    └── pal_neog.tre   <- copied into pal_neog/\n\n"
        "In this example, 'pal_neog.tre' is used as the tree for all FASTAs in 'pal_neog/'.\n\n" \
        "  • Use --treedir (and --tree-ext if necessary) to use a tree for each FASTA by matching basename from a directory.\n" \
        "       example: for FASTA OG0001.fa and --tree-dir trees/ --tree-ext .tre, we look for trees/OG0001.tre.\n"
        "  • When --tree-dir-search is set, and run-dirs is processing subdirectory X, the script searches\n" 
        "    for trees under PATH(--tree-dir)/X/ and matches by FASTA basename (e.g., OG0001_aln.fas → OG0001_aln.fa.tre.rtd).\n" 
        "    If the directory PATH(--tree-dir)/X/ is not found, that subdirectory is skipped."
    ),
    formatter_class=argparse.RawDescriptionHelpFormatter
)
    p3.add_argument('-e', '--fa-ext', default='.fa', help='FASTA extension (default: .fa)')
    p3.add_argument('-r', '--resume', action='store_true', help='Skip if .FigTree.tre exists in subdir')
    p3.add_argument('--ctl-template', help='Override .ctl template')
    p3.add_argument('--wag-file', help='Override WAG matrix file')
    p3.add_argument('--tree-dir', help='Directory containing per-gene trees; match by FASTA basename')
    p3.add_argument('--tree-ext', default='.tre', help='Tree filename extension in --tree-dir (default: .tre)')
    p3.add_argument(
        '--sanitize',
        nargs='?', const='',
        help=("Sanitize per-gene trees from --tree-dir using letters "
              "(l=keep-lengths, s=keep-support, p=allow-polytomies, h=add-header, c=no-collapse). "
              "Example: '--sanitize sc'")
    )
    p3.add_argument(
    "--tree-dir-search",
    action="store_true",
    help=("Search for trees under SUBTREE = PATH(--tree-dir)/<subdir_name>/ "
          "when matching trees by FASTA basename. If the matching directory "
          "does not exist, the subdirectory is skipped. Requires --tree-dir.")
)
    p3.set_defaults(func=cmd_run_dirs)

