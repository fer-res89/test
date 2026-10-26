#!/usr/bin/env python3
"""
aeon — one-stop CLI for:
  • bucket-files        Bucket gene trees by quartet support and copy matching FASTAs
  • process             Run MCMCTree + CodeML on a single FASTA
  • pipe-all            Run the pipeline on all FASTAs in cwd
  • run-dirs            Run the pipeline across subdirectories
  • extract-dates       Extract 95% HPD ages from NEXUS .tre files
  • rearrange-triplet   Produce (g1,g3)|g2 and (g2,g3)|g1 Newicks preserving original root split/order
  • sanitize-tree       Sanitize a Newick file (collapse unary, drop numeric internal labels, etc.)

Any subcommand can be run once per subdirectory with the global --per-dir
flag (see `aeon --help`).

Dependencies:
  - dendropy, Biopython (Bio.Phylo), ete3
  - External tools in $PATH: mcmctree, codeml

Adding a new subcommand: create a new module in this package exposing a
`register(sub)` function that adds its subparser (see any existing module
for the pattern), then add it to SUBCOMMAND_MODULES below.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from . import bucketing, process, extract_dates, rearrange_triplet, sanitize_tree
except ModuleNotFoundError as e:
    sys.exit(
        f"aeon is missing a required dependency ({e.name}).\n"
        f"Install it with: pip install 'phylotools[aeon]'"
    )

# Modules registered as `aeon <name>` subcommands, in help-listing order.
SUBCOMMAND_MODULES = (
    bucketing,
    process,
    extract_dates,
    rearrange_triplet,
    sanitize_tree,
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Unified phylogenetics toolkit (bucketing, pipeline, HPD extraction, rearrangement)",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # -------- Global per-directory mode --------
    p.add_argument(
        "--per-dir",
        action="store_true",
        help=(
            "If set, run the selected COMMAND once inside each immediate subdirectory of "
            "--per-dir-root. For each run, the working directory is that subdirectory, so "
            "relative paths in other options are interpreted relative to it."
        ),
    )
    p.add_argument(
        "--per-dir-root",
        default=".",
        help=(
            "Root directory whose immediate subdirectories are processed when --per-dir is "
            "given (default: current directory)."
        ),
    )
    p.add_argument(
        "--per-dir-include-hidden",
        action="store_true",
        help="Include subdirectories whose names start with '.' when using --per-dir.",
    )
    p.add_argument(
        "--per-dir-jobs",
        type=int,
        default=1,
        help=(
            "When using --per-dir, maximum number of subdirectories to process in parallel "
            "(default 1 = sequential)."
        ),
    )
    p.add_argument(
        "--per-dir-recursive-target",
        metavar="NAME",
        help=(
            "If set, instead of running in each immediate subdirectory of --per-dir-root, "
            "search *inside* each of those subdirectories for subdirectories whose name is "
            "exactly NAME, and run the selected COMMAND in each such match.\n\n"
            "Example: with --per-dir-root . and --per-dir-recursive-target runs,\n"
            "aeon will run the COMMAND inside every directory named 'runs' that "
            "is found anywhere under each immediate child of --per-dir-root."
        ),
    )

    sub = p.add_subparsers(dest="cmd", metavar="COMMAND", required=True)

    for mod in SUBCOMMAND_MODULES:
        mod.register(sub)

    return p


# ============================================================================
# --per-dir execution
# ============================================================================

def _per_dir_worker(func, args_dict, subdir_str):
    """
    Worker for per-dir parallelization.

    - func: the subcommand function to run (e.g. cmd_bucket_files)
    - args_dict: dict(args) from the top-level parser
    - subdir_str: directory in which to run the command
    """
    # Each worker is a separate process; change cwd locally
    os.chdir(subdir_str)

    # Build a fresh Namespace for this run
    inner_dict = args_dict.copy()
    # Avoid recursion inside workers
    inner_dict["per_dir"] = False
    # In case user refers to it, make root "." inside the subdir
    inner_dict["per_dir_root"] = "."

    inner_args = argparse.Namespace(**inner_dict)

    # Run the actual command
    func(inner_args)
    return subdir_str


def run_per_dir(args):
    """
    If --per-dir was requested, run args.func(args) once per target directory.

    Targets are:

    - By default: each immediate subdirectory of --per-dir-root
      (optionally skipping hidden dirs).

    - If --per-dir-recursive-target NAME is given:
      For each immediate subdirectory under --per-dir-root, search recursively
      for subdirectories whose *basename* is NAME, and use those as targets.
    """
    root = Path(args.per_dir_root).resolve()
    if not root.is_dir():
        sys.exit(f"[ERROR] --per-dir-root is not a directory: {root}")

    include_hidden = getattr(args, "per_dir_include_hidden", False)
    recursive_target = getattr(args, "per_dir_recursive_target", None)

    # Immediate children under root (these define the "outer" level)
    outer_dirs = sorted(
        d for d in root.iterdir()
        if d.is_dir() and (include_hidden or not d.name.startswith("."))
    )

    if not outer_dirs:
        print(f"[INFO] No subdirectories found under {root}")
        return

    target_dirs: list[Path] = []

    if recursive_target:
        # CONDITIONAL RECURSIVE MODE:
        # For each outer dir, find all subdirs named exactly recursive_target.
        for outer in outer_dirs:
            # NOTE: We allow matches even if they live under hidden intermediate
            # dirs; only the outer-level inclusion is controlled by include_hidden.
            for sub in outer.rglob(recursive_target):
                if sub.is_dir():
                    target_dirs.append(sub)
        if not target_dirs:
            print(
                f"[INFO] --per-dir-recursive-target={recursive_target!r} was given, "
                f"but no matching subdirectories were found under {root}"
            )
            return
    else:
        # FLAT MODE: just run in each immediate child
        target_dirs = outer_dirs

    # How many workers?
    max_workers = getattr(args, "per_dir_jobs", 1) or 1
    max_workers = max(1, min(max_workers, len(target_dirs)))

    func = args.func
    args_dict = vars(args).copy()

    if max_workers == 1 or len(target_dirs) == 1:
        # Sequential mode
        for sub in target_dirs:
            print(f"[per-dir] Running '{args.cmd}' in {sub}")
            _per_dir_worker(func, args_dict, str(sub))
    else:
        print(
            f"[per-dir] Running '{args.cmd}' in {len(target_dirs)} directories "
            f"with up to {max_workers} parallel workers."
        )
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            future_to_sub = {
                ex.submit(_per_dir_worker, func, args_dict, str(sub)): sub
                for sub in target_dirs
            }
            for fut in as_completed(future_to_sub):
                sub = future_to_sub[fut]
                try:
                    fut.result()
                    print(f"[per-dir] Completed in {sub}")
                except Exception as e:
                    # Report failure but don't stop others
                    print(f"[per-dir] ERROR in {sub}: {e}", file=sys.stderr)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if getattr(args, "per_dir", False):
        run_per_dir(args)
    else:
        rc = args.func(args)
        if rc is False:          # subcommands that report failure (process, pipe-all, run-dirs)
            sys.exit(1)


if __name__ == "__main__":
    main()
