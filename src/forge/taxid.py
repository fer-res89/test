"""The `forge get-taxid` subcommand: build a taxid_map.tsv from GPFF
files.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from typing import List, Optional, Tuple

try:
    from Bio import SeqIO  # type: ignore
except Exception:
    SeqIO = None  # type: ignore

def expand_inputs(patterns: List[str]) -> List[str]:
    files: List[str] = []
    for pat in patterns:
        matches = glob.glob(pat)
        if matches:
            files.extend(matches)
        elif os.path.isfile(pat):
            files.append(pat)
        else:
            print(f"Warning: '{pat}' did not match any files", file=sys.stderr)
    return sorted(set(files))

def extract_taxid_map(gpff_files: List[str], output_path: str) -> None:
    if SeqIO is None:
        sys.exit("BioPython is required for 'get-taxid'. Please install biopython.")
    entries: List[Tuple[str, str]] = []
    for gpff in gpff_files:
        for rec in SeqIO.parse(gpff, 'genbank'):
            taxid: Optional[str] = None
            for feat in rec.features:
                if feat.type == 'source':
                    for xref in feat.qualifiers.get('db_xref', []):
                        if xref.startswith('taxon:'):
                            taxid = xref.split(':', 1)[1]
                            break
                if taxid:
                    break
            if taxid:
                entries.append((rec.id, taxid))
            else:
                print(f"Warning: no taxon ID for record '{rec.id}' in {gpff}", file=sys.stderr)
    with open(output_path, 'w') as out_h:
        for seqid, taxid in entries:
            out_h.write(f"{seqid}\t{taxid}\n")
    print(f"Wrote {len(entries)} entries to {output_path}", file=sys.stderr)


def register(sub, formatter_class) -> None:
    """Register the `get-taxid` subcommand."""
    def run_get_taxid(args: argparse.Namespace) -> None:
        gpff_files = expand_inputs(args.inputs)
        if not gpff_files:
            sys.exit("Error: No GPFF files found.")
        extract_taxid_map(gpff_files, args.output)

    sp_tax = sub.add_parser(
        "get-taxid", 
        help=("Build taxid_map.tsv from GPFF files.\n"
              "(run 'forge get-taxid -h' for options)"
        ),
    formatter_class=formatter_class,
    )
    sp_tax.add_argument("-i", "--inputs", nargs="+", required=True, help="GPFF files or globs (e.g., '*.gbff')")
    sp_tax.add_argument("-o", "--output", default="taxid_map.tsv", help="Output TSV (default: taxid_map.tsv)")
    sp_tax.set_defaults(func=run_get_taxid)
