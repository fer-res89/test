"""The `forge dedupe` subcommand and the FASTA de-duplication-by-taxon
helpers reused by the `filter` subcommand.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set

try:
    from Bio import SeqIO  # type: ignore
except Exception:
    SeqIO = None  # type: ignore

from phylotools_common.fileutils import ensure_dir
from ._common import FASTA_GLOBS, species_from_species_tree

def _norm(s: str) -> str:
    s = s.lower()
    s = re.sub(r'[^a-z0-9]+', '_', s)
    s = re.sub(r'_+', '_', s).strip('_')
    return s

def build_species_norms(species_list: List[str]) -> Dict[str, List[str]]:
    spec_norms: Dict[str, List[str]] = {}
    for sp in species_list:
        variants = {
            sp,
            sp.replace(' ', '_'),
            sp.replace('_', ' '),
            sp.replace('-', '_'),
            sp.replace('_', '-'),
            sp.replace('-', ' '),
        }
        parts = re.split(r'[ _-]+', sp.strip())
        if len(parts) >= 3:
            two = f"{parts[0]}_{parts[1]}"
            variants |= {two, two.replace('_', ' '), two.replace('_', '-')}
        normed = sorted({_norm(v) for v in variants if v.strip()}, key=len, reverse=True)
        if _norm(sp) not in normed:
            normed.append(_norm(sp))
        spec_norms[sp] = normed
    return spec_norms

def choose_species_from_header_partial(header: str,
                                       species_norms: Dict[str, List[str]]) -> Optional[str]:
    h = _norm(header)
    best_sp, best_len = None, -1
    for canon, aliases in species_norms.items():
        for a in aliases:
            if a and a in h and len(a) > best_len:
                best_sp, best_len = canon, len(a)
    return best_sp

def fasta_paths(directory: Path) -> List[Path]:
    paths: List[Path] = []
    for pat in FASTA_GLOBS:
        paths.extend(sorted(directory.glob(pat)))
    return [p for p in paths if "_clean" not in p.stem and "_edit" not in p.stem]

def dedupe_file_inplace(fasta_path: Path,
                        species_norms: Dict[str, List[str]]) -> int:
    if SeqIO is None:
        sys.exit("BioPython is required for FASTA operations. Install biopython.")
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    seen: Set[str] = set()
    kept = []
    for rec in records:
        header = rec.id or ""
        sp = choose_species_from_header_partial(header, species_norms)
        key = sp if sp is not None else f"__NO_SPECIES__::{_norm(header)}"
        if key in seen:
            continue
        seen.add(key)
        kept.append(rec)
    tmp = fasta_path.with_suffix(fasta_path.suffix + ".tmp")
    with tmp.open("w") as fh:
        SeqIO.write(kept, fh, "fasta")
    tmp.replace(fasta_path)
    return len(kept)

def run_dedupe_dir_inplace(fasta_dir: Path,
                           taxa_list: Optional[List[str]],
                           species_tree_list: List[str]) -> int:
    species = taxa_list if taxa_list else species_tree_list
    if not species:
        sys.exit("No species list available: provide --taxa-file or ensure Species_Tree/SpeciesTree_rooted.txt exists.")
    species_norms = build_species_norms(species)
    n_files = 0
    for p in fasta_paths(fasta_dir):
        dedupe_file_inplace(p, species_norms)
        n_files += 1
    return n_files


def register(sub, formatter_class) -> None:
    """Register the `dedupe` subcommand."""
    # dedupe
    def run_dedupe_cmd(args: argparse.Namespace) -> None:
        base = Path(args.base_dir)

        # Dedupe can operate on any FASTA directory; don't require OrthoFinder layout here.
        fasta_dir = base / (args.fasta_dir if args.fasta_dir else Path("Orthogroup_Sequences/all_taxa_fas"))
        ensure_dir(fasta_dir)

        taxa_list = None
        if args.taxa_file and os.path.isfile(args.taxa_file):
            with open(args.taxa_file) as f:
                taxa_list = [line.strip() for line in f if line.strip()]

        # Fallback: try to read species from a species tree file if it exists.
        species_fallback: List[str] = []
        cand_paths = [
            base / "Species_Tree" / "SpeciesTree_rooted.txt",
            base / "Species_Tree" / "SpeciesTree.txt",
            base / "SpeciesTree_rooted.txt",
            base / "SpeciesTree.txt",
        ]
        for p in cand_paths:
            if p.is_file():
                species_fallback = species_from_species_tree(p)
                break

        n = run_dedupe_dir_inplace(fasta_dir, taxa_list, species_fallback)
        print(f"Deduped {n} file(s) in: {fasta_dir}")

    sp_ded = sub.add_parser(
        "dedupe", 
        help=("Deduplicate FASTAs in-place by taxon (first match kept).\n"
              "(run 'forge dedupe -h' for options)"
        ),
    formatter_class=formatter_class,
    ) 
    sp_ded.add_argument("--base-dir", type=Path, default=Path("."),
                        help="Directory wiht Orthofinder results. Orthobrew will search for Orthogroup_Sequences/all_taxa_fas and Species_Tree/SpeciesTree_rooted.txt if no other arguments are used")
    sp_ded.add_argument("--fasta-dir", type=Path, default=None,
                        help="Target FASTA directory (default: Orthogroup_Sequences/all_taxa_fas)")
    sp_ded.add_argument("--taxa-file", type=str, default=None, help="Taxa list for dedupe (one per line).")
    sp_ded.set_defaults(func=run_dedupe_cmd)
