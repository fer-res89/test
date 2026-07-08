"""The `extract-dates` subcommand: extract 95% HPD ages from NEXUS .tre
files produced by MCMCTree.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

from Bio import Phylo

def _ed_load_quartets_tsv(tsv_path: str) -> dict[str, dict[str, set[str]]]:
    """
    Load quartets.tsv into a dict: name -> { "A": set(...), "B": set(...), "C": set(...), "D": set(...) }.
    Same format as bucket-files.
    """
    import csv
    qmap: dict[str, dict[str, set[str]]] = {}
    with open(tsv_path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            nm = (row.get("name") or "").strip()
            if not nm:
                continue

            def s(field: str) -> set[str]:
                raw = (row.get(field) or "")
                return {t.strip() for t in raw.split(",") if t.strip()}

            qmap[nm] = {
                "A": s("A_set"),
                "B": s("B_set"),
                "C": s("C_set"),
                "D": s("D_set"),
            }
    return qmap


def _ed_build_in_ext_from_quartet(qrec: dict[str, set[str]]) -> tuple[list[str], list[str]]:
    """
    Given a quartet record with keys A,B,C,D (each a set of taxa),
    return (tip_in, tip_ex) as *sorted* lists:

      tip_in = A ∪ B
      tip_ex = A ∪ B ∪ C

    This matches your desired behavior:
      - A_set + B_set : -tp -in
      - A_set + B_set + C_set : -tp -ex
    """
    A, B, C = qrec["A"], qrec["B"], qrec["C"]
    tip_in = sorted(A | B)
    tip_ex = sorted(A | B | C)
    return tip_in, tip_ex



def _ed_run_on_dir_simple(
    data_dir: Path,
    out_dir: Path,
    tips: list[str],
    args,
):
    """
    Original extract-dates behavior: ONE CSV per data_dir, ONE tip-set.

    - data_dir: where to look for tree files
    - out_dir: where to write the CSV
    - tips: target tip set (from --target-tips or similar)
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    prefix = (
        'dates_int_'
        if args.mode == 'int' else
        'dates_ext_'
        if args.mode == 'ext' else
        'dates_'
    )

    out_csv = out_dir / f"{prefix}{data_dir.name}.csv"

    files = [f for f in os.listdir(data_dir) if f.endswith(args.extension)]
    skipped_tree, skipped_hpd = [], []

    with open(out_csv, 'w') as out:
        for fname in files:
            fpath = data_dir / fname
            ages = _get_ages(str(fpath), tips)
            if ages is None:
                try:
                    Phylo.read(str(fpath), 'nexus')
                except ValueError:
                    skipped_tree.append(fname)
                else:
                    skipped_hpd.append(fname)
            else:
                lo, hi, avg = ages
                base = fname.rsplit(args.extension, 1)[0]
                out.write(f"{base}, {lo:.6f}, {hi:.6f}, {avg:.6f}\n")

    print(f"Wrote results to {out_csv}")
    if skipped_tree:
        print("Skipped (no trees):", ', '.join(skipped_tree))
    if skipped_hpd:
        print("Skipped (no HPD):", ', '.join(skipped_hpd))




def _ed_run_on_dir_quartet(
    data_dir: Path,
    out_dir: Path,
    tip_in: list[str],
    tip_ext: list[str],
    args,
):
    """
    Quartet-aware mode:
      - write dates_int_<dir>.csv using tip_in  (A∪B)
      - write dates_ext_<dir>.csv using tip_ext (A∪B∪C)
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    files = [f for f in os.listdir(data_dir) if f.endswith(args.extension)]

    # -------- INT: A∪B --------
    out_int = out_dir / f"dates_int_{data_dir.name}.csv"
    skipped_tree_int, skipped_hpd_int = [], []

    with open(out_int, 'w') as out:
        for fname in files:
            fpath = data_dir / fname
            ages = _get_ages(str(fpath), tip_in)
            if ages is None:
                try:
                    Phylo.read(str(fpath), 'nexus')
                except ValueError:
                    skipped_tree_int.append(fname)
                else:
                    skipped_hpd_int.append(fname)
            else:
                lo, hi, avg = ages
                base = fname.rsplit(args.extension, 1)[0]
                out.write(f"{base}, {lo:.6f}, {hi:.6f}, {avg:.6f}\n")

    print(f"[INT] Wrote results to {out_int}")
    if skipped_tree_int:
        print("[INT] Skipped (no trees):", ', '.join(skipped_tree_int))
    if skipped_hpd_int:
        print("[INT] Skipped (no HPD):", ', '.join(skipped_hpd_int))

    # -------- EXT: A∪B∪C --------
    out_ext = out_dir / f"dates_ext_{data_dir.name}.csv"
    skipped_tree_ext, skipped_hpd_ext = [], []

    with open(out_ext, 'w') as out:
        for fname in files:
            fpath = data_dir / fname
            ages = _get_ages(str(fpath), tip_ext)
            if ages is None:
                try:
                    Phylo.read(str(fpath), 'nexus')
                except ValueError:
                    skipped_tree_ext.append(fname)
                else:
                    skipped_hpd_ext.append(fname)
            else:
                lo, hi, avg = ages
                base = fname.rsplit(args.extension, 1)[0]
                out.write(f"{base}, {lo:.6f}, {hi:.6f}, {avg:.6f}\n")

    print(f"[EXT] Wrote results to {out_ext}")
    if skipped_tree_ext:
        print("[EXT] Skipped (no trees):", ', '.join(skipped_tree_ext))
    if skipped_hpd_ext:
        print("[EXT] Skipped (no HPD):", ', '.join(skipped_hpd_ext))



def _ed_build_pair_in_ext_from_quartet(qrec: dict, which: str):
    """
    Given a quartet record qrec with keys 'A','B','C','D',
    and a topology label like 'AB|CD', 'AC|BD', or 'AD|BC',
    return (tip_in, tip_ext), where:

      - tip_in is the ingroup pair:
          AB|CD → A∪B
          AC|BD → A∪C
          AD|BC → B∪C

      - tip_ext is always A∪B∪C (deeper comparison node).
    """
    A = set(qrec["A"])
    B = set(qrec["B"])
    C = set(qrec["C"])
    # D is not needed for these ages, but keep for clarity
    _D = set(qrec.get("D", []))

    if which == "AB|CD":
        ingroup = A | B
    elif which == "AC|BD":
        ingroup = A | C
    elif which == "AD|BC":
        ingroup = B | C
    else:
        # Fallback: default to AB|CD behaviour if something weird happens
        ingroup = A | B

    ext = A | B | C
    return ingroup, ext


# =========================
# bucket-files
# =========================


def _get_ages(tree_file, target_tips):
    try:
        tree = Phylo.read(tree_file, 'nexus')
    except ValueError:
        return None
    mrca = tree.common_ancestor(target_tips)
    if not mrca or not mrca.comment:
        return None
    parts = mrca.comment.strip('[]').split('&95%HPD=')
    if len(parts) != 2:
        return None
    lo, hi = map(float, parts[1].strip('{}').split(','))
    return lo, hi, (lo + hi) / 2

def cmd_extract_dates(args):
    """
    Extended extract-dates:

    - Default (no quartets / no auto / no batch):
        * exactly original behavior: one CSV per directory using --target-tips and --mode

    - With -qt/--quartets-tsv + -auto/--auto-quartets (non-batch):
        * interpret cwd.name as quartet name
        * infer A,B,C from quartets.tsv
        * write:
            dates_int_<cwd>.csv  (A∪B)
            dates_ext_<cwd>.csv  (A∪B∪C)

    - With -qt/--quartets-tsv + -b/--batch + -out/--out-dir:
        * from CWD, for each quartet 'name' that has a matching directory:
            data_dir = ./<name>
            out_dir  = <out-dir>/<name>/
          write:
            dates_int_<name>.csv
            dates_ext_<name>.csv
    """

    # ---- Enforce target-tips requirement only when needed ----
    if not getattr(args, "quartets_tsv", None):
        # No quartets file → must provide target-tips
        if not getattr(args, "target_tips", None):
            sys.exit(
                "[extract-dates] Error: --target-tips is required when no --quartets-tsv is given."
            )
    else:
        # quartets TSV exists: only require target-tips if NOT using auto/batch
        if not (getattr(args, "auto_quartets", False) or getattr(args, "batch", False)):
            if not getattr(args, "target_tips", None):
                sys.exit(
                    "[extract-dates] Error: --target-tips is required unless using "
                    "--auto-quartets or --batch with --quartets-tsv."
                )

    # ---- Load quartets map if provided ----
    qmap: dict[str, dict[str, set[str]]] = {}
    if getattr(args, "quartets_tsv", None):
        qmap = _ed_load_quartets_tsv(args.quartets_tsv)

    # ---- Output root ----
    out_root = Path(args.out_dir) if getattr(args, "out_dir", None) else None

    # =====================
    # MODE A: non-batch
    # =====================
    if not getattr(args, "batch", False):
        data_dir = Path.cwd()

        # Base tips from CLI (legacy/simple behavior)
        if getattr(args, "target_tips", None):
            tips = [t.strip() for t in args.target_tips.split(",") if t.strip()]
        else:
            tips = []

        # No auto-quartets or no quartets TSV → original behavior
        if not (getattr(args, "auto_quartets", False) and qmap):
            # choose out_dir (out_root if given, else current dir)
            out_dir = out_root if out_root is not None else data_dir
            _ed_run_on_dir_simple(data_dir, out_dir, tips, args)
            return

        # auto-quartets mode with quartets.tsv
        qname = data_dir.name
        qrec = qmap.get(qname)
        if qrec is None:
            print(
                f"[extract-dates] --auto-quartets: no quartet named '{qname}' "
                f"in {args.quartets_tsv}; falling back to simple mode.",
                file=sys.stderr,
            )
            out_dir = out_root if out_root is not None else data_dir
            _ed_run_on_dir_simple(data_dir, out_dir, tips, args)
            return

        tip_in, tip_ext = _ed_build_in_ext_from_quartet(qrec)
        #print(
        #    f"[extract-dates] Quartet '{qname}': "
        #    f"{len(tip_in)} internal tips (A∪B), {len(tip_ext)} external tips (A∪B∪C).",
        #    file=sys.stderr,
        #)

        # out_dir: <out-root>/<quartet>/ if out_root is given, else current dir
        out_dir = out_root / qname if out_root is not None else data_dir
        _ed_run_on_dir_quartet(data_dir, out_dir, tip_in, tip_ext, args)
        return

    # =====================
    # MODE B: batch over quartet directories
    # =====================
    if not getattr(args, "batch", False):
        # fall back to your non-batch logic (what you already have)
        # e.g. old cmd_extract_dates body
        prefix = 'dates_int_' if args.mode == 'int' else 'dates_ext_' if args.mode == 'ext' else 'dates_'
        cwd = Path.cwd()
        tips = [t.strip() for t in args.target_tips.split(',')]
        files = [f for f in os.listdir(cwd) if f.endswith(args.extension)]
        out_csv = f"{prefix}{cwd.name}.csv"
        skipped_tree, skipped_hpd = [], []
        with open(out_csv, 'w') as out:
            for f in files:
                ages = _get_ages(f, tips)
                if ages is None:
                    try:
                        Phylo.read(f, 'nexus')
                    except ValueError:
                        skipped_tree.append(f)
                    else:
                        skipped_hpd.append(f)
                else:
                    lo, hi, avg = ages
                    base = f.rsplit(args.extension, 1)[0]
                    out.write(f"{base}, {lo:.6f}, {hi:.6f}, {avg:.6f}\n")
        print(f"Wrote results to {out_csv}")
        if skipped_tree:
            print("Skipped (no trees):", ', '.join(skipped_tree))
        if skipped_hpd:
            print("Skipped (no HPD):", ', '.join(skipped_hpd))
        return

    # ---- from here downward: true batch mode ----

    if not args.quartets_tsv:
        sys.exit("[extract-dates] --batch requires --quartets-tsv to be set.")

    if not args.out_dir:
        sys.exit("[extract-dates] --batch requires --out-dir to specify where CSVs go.")

    root = Path.cwd()
    out_root = Path(args.out_dir).resolve()

    # Load all quartet records from the TSV into a dict like {name: qrec}
    
    quartets = _qb_load_quartets(None, args.quartets_tsv)  # same helper as bucket-files

    for qrec in quartets:
        qname = qrec["name"]
        qroot = root / qname
        if not qroot.is_dir():
            continue

        # Build topology→dirname mapping consistent with bucket-files
        if "_qb_bucket_dir_map" in globals():
            dir_map = _qb_bucket_dir_map(qrec, getattr(args, "focal", None))
        else:
            dir_map = {
                "AB|CD": "AB_CD",
                "AC|BD": "AC_BD",
                "AD|BC": "AD_BC",
                "unresolved": "unresolved",
                "insufficient": "insufficient",
            }


        print(f"[DEBUG] Processing quartet: {qname}", file=sys.stderr)
        print(f"[DEBUG] dir_map = {dir_map}", file=sys.stderr)
        print(
            f"[DEBUG] subdirectories = "
            f"{[d.name for d in qroot.iterdir() if d.is_dir()]}",
            file=sys.stderr,
        )            


        # Invert it: dirname → topology (only the informative three)
        topo_for_dir = {}
        for topo in ("AB|CD", "AC|BD", "AD|BC"):
            dname = dir_map.get(topo)
            if dname:
                topo_for_dir[dname] = topo

        # Subdirectories of the quartet dir are the data dirs
        subdirs = sorted(d for d in qroot.iterdir() if d.is_dir())
        if not subdirs:
            continue

        q_out_root = out_root / qname

        # Iterate each data dir (e.g. archo_lep, turt_lep, arche_lep)
        for data_dir in subdirs:
            topo = topo_for_dir.get(data_dir.name)
            if topo is None:
                # e.g. unresolved/insufficient or some extra dir – skip
                print(
                    f"[extract-dates] skipping {data_dir} (no topology mapping for this dir name)",
                    file=sys.stderr,
                )
                continue

            tip_in, tip_ext = _ed_build_pair_in_ext_from_quartet(qrec, topo)
            print(
                f"[extract-dates] Quartet '{qname}', dir '{data_dir.name}': "
                f"topology {topo}, |in|={len(tip_in)}, |ext|={len(tip_ext)}",
                file=sys.stderr,
            )

            _ed_run_on_dir_quartet(
                data_dir=data_dir,
                out_dir=q_out_root,
                tip_in=tip_in,
                tip_ext=tip_ext,
                args=args,
            )




# =========================
# Rearrangement
# =========================



def register(sub) -> None:
    """Register the `extract-dates` subcommand."""
    # extract-dates
    p4 = sub.add_parser('extract-dates', help='Extract 95%% HPD ages from NEXUS .tre files')
    p4.add_argument('-in', dest='mode', action='store_const', const='int', help='Use dates_int_ prefix')
    p4.add_argument('-ex', dest='mode', action='store_const', const='ext', help='Use dates_ext_ prefix')
    p4.add_argument('-tp', '--target-tips', required=False, help='Comma-separated tip labels')
    p4.add_argument('-e', '--extension', default='.tre', help='Tree file extension (default: .tre)')
    p4.add_argument(
        "-qt", "--quartets-tsv",
        help="TSV with columns: name, A_set, B_set, C_set, D_set (same format as bucket-files). "
             "If set together with --auto-quartets, A_set∪B_set is used as ingroup, "
             "A_set∪B_set∪C_set as exclude."
    )
    p4.add_argument(
        "-auto", "--auto-quartets",
        action="store_true",
        help="Look for quartet definitions in --quartets-tsv whose 'name' matches the "
             "quartet directory name. For each quartet dir, use A_set∪B_set as for -in "
             "and A_set∪B_set∪C_set for -ex when extracting dates."
    )
    p4.add_argument(
        "-out", "--out-dir",
        help="Optional root directory for outputs. If provided, CSVs for data directory DIR "
             "under quartet Q will be written to <out-dir>/Q/DIR/. If omitted, behave "
             "like the current extract-dates behavior."
    )
    p4.add_argument(
        "-b", "--batch",
        action="store_true",
        help="Batch mode: starting from CWD, look for directories whose name matches a "
             "quartet 'name' in --quartets-tsv, and inside each such directory process "
             "its immediate subdirectories as data dirs."
    )
    p4.set_defaults(func=cmd_extract_dates)

