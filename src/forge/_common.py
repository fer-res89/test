"""Constants and tree/file helpers shared across forge subcommands.

Uses Bio.Phylo for newick parsing (kept as-is per project decision; not
consolidated with aeon's dendropy-based tree handling).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import List, Set

try:
    from Bio import Phylo  # type: ignore
except Exception:
    Phylo = None  # type: ignore

from phylotools_common.fileutils import ensure_dir, sanitize_id

TREE_FILE_RE = re.compile(r"^OG\d+_tree.*\.txt$", re.IGNORECASE)
TREE_OG_RE = re.compile(r"^(?P<og>.+?)_tree.*\.txt$", re.IGNORECASE)
FASTA_GLOBS = ("*.fa", "*.fas", "*.fasta", "*.FA", "*.FAS", "*.FASTA")
LABEL_NEWICK_RE = re.compile(r"^\s*(?P<label>[^:]+?)\s*:\s*(?P<newick>\(.*)\s*$")


def count_tree_files(d: Path) -> int:
    return sum(1 for p in d.iterdir() if p.is_file() and TREE_FILE_RE.match(p.name))


def count_fasta_files(d: Path) -> int:
    n = 0
    for pat in FASTA_GLOBS:
        n += len(list(d.glob(pat)))
    return n


def materialize_trees_if_needed(
    trees_path: Path,
    work_dir: Path,
    on_progress=None,
    wanted_ids=None,
) -> Path:
    """
    If trees_path is a directory: return it.
    If it's a file with lines like 'LABEL: <newick>': split into per-tree files:
        <LABEL>_tree.txt
    LABEL can be OG..., tree1, gene_tree1, etc.

    on_progress: optional callable invoked as on_progress(n_written) as trees
    are written, so a caller can display live progress during a large split.

    wanted_ids: optional Set[str] of raw labels to keep. When given, any line
    whose LABEL is not in this set is skipped entirely (not written to disk).
    This matters because callers (e.g. forge filter) typically only want the
    all-taxa-orthogroup subset -- writing every orthogroup in a big
    Resolved_Gene_Trees.txt to disk just to discard most of them downstream
    is pure waste for large runs. Matching is done on the RAW label (before
    sanitize_id()), since that's what identifiers like og_ids from
    Orthogroups.GeneCount.tsv look like; sanitize_id() is a no-op for
    ordinary OG-style names anyway.
    """
    if trees_path.is_dir():
        return trees_path

    ensure_dir(work_dir)

    n_written = 0
    n_scanned = 0
    used = set()

    with trees_path.open("r") as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue

            m = LABEL_NEWICK_RE.match(line)
            if not m:
                raise ValueError(
                    f"Line {line_no} doesn't look like 'LABEL: <newick>': {line[:160]}"
                )

            label_raw = m.group("label").strip()
            n_scanned += 1

            if wanted_ids is not None and label_raw not in wanted_ids:
                continue

            newick = m.group("newick").strip()
            if not newick.endswith(";"):
                newick += ";"

            label = sanitize_id(label_raw)

            # Avoid collisions after sanitizing (rare but possible)
            out_id = label
            k = 2
            while out_id in used:
                out_id = f"{label}__{k}"
                k += 1
            used.add(out_id)

            out = work_dir / f"{out_id}_tree.txt"
            out.write_text(newick + "\n")
            n_written += 1
            if on_progress is not None and (n_written % 50 == 0):
                on_progress(n_written)

    if on_progress is not None:
        on_progress(n_written)

    if n_written == 0:
        if wanted_ids is not None and n_scanned > 0:
            raise ValueError(
                f"None of the {n_scanned} tree(s) in {trees_path} matched any "
                f"of the {len(wanted_ids)} requested id(s). Check that labels "
                f"in the trees file (e.g. 'OG0000001: (...);') match the "
                f"orthogroup ids you're filtering to exactly."
            )
        raise ValueError(f"No trees written from {trees_path}")

    return work_dir


def species_from_species_tree(tree_path: Path) -> List[str]:
    if not tree_path.is_file():
        return []
    if Phylo is None:
        sys.exit("BioPython is required for species tree parsing. Install biopython.")
    tree = Phylo.read(str(tree_path), "newick")
    seen: Set[str] = set()
    species: List[str] = []
    for term in tree.get_terminals():
        nm = (term.name or "").strip()
        if nm and nm not in seen:
            seen.add(nm)
            species.append(nm)
    return species
