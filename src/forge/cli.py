#!/usr/bin/env python3
"""
forge — filter/QC OrthoFinder output, clean FASTA/newick labels,
deduplicate FASTAs by taxon, append best-BLAST-hit sequences, and build
taxid maps from GPFF files.

Examples:
  forge filter --base-dir .
  forge filter --base-dir . --dedupe
  forge filter --base-dir . --clean-labels taxa.txt
  forge filter --base-dir . --dedupe --clean-labels taxa.txt --inplace

  forge clean-labels -i 'Gene_Trees/all_taxa_orthogroups/*.txt' -f newick -t taxa.txt --inplace
  forge blast-species -e .fasta -db target_db -tl blastn -t 8
  forge get-taxid -i '*.gbff' -o taxid_map.tsv

Adding a new subcommand: create a new module in this package exposing a
`register(sub, formatter_class)` function that adds its subparser (see any
existing module for the pattern), then add it to SUBCOMMAND_MODULES below.
"""
from __future__ import annotations

import argparse
import sys

try:
    from . import filter_pipeline, clean_labels, dedupe, blast_species, miniprot_recover, taxid
except ModuleNotFoundError as e:
    sys.exit(
        f"forge is missing a required dependency ({e.name}).\n"
        f"Install it with: pip install 'phylotools[forge]'"
    )

# Modules registered as `forge <name>` subcommands, in help-listing order.
SUBCOMMAND_MODULES = (
    filter_pipeline,
    clean_labels,
    dedupe,
    blast_species,
    miniprot_recover,
    taxid,
)


class HelpFmt(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    """Show defaults and keep newlines in help/epilog text."""
    pass


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Unified phylogenomic toolkit (ortho-groups filter, label "
            "cleaning, BLAST per-species, GPFF taxid map)"
        ),
        formatter_class=HelpFmt,
        epilog=__doc__,
    )
    sub = p.add_subparsers(dest="tool", required=True)

    for mod in SUBCOMMAND_MODULES:
        mod.register(sub, HelpFmt)

    return p


def main(argv=None) -> int:
    """Console entry point for the forge CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    try:
        args.func(args)
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
