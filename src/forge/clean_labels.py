"""The `forge clean-labels` subcommand and the FASTA/newick label
cleaning helpers reused by the `filter` subcommand.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from typing import List, Optional, Tuple

try:
    from Bio import SeqIO, Phylo  # type: ignore
except Exception:
    SeqIO = None  # type: ignore
    Phylo = None  # type: ignore

def load_taxon_filter(taxon_arg: str) -> Tuple[Optional[List[str]], Optional[re.Pattern]]:
    if os.path.isfile(taxon_arg):
        with open(taxon_arg) as f:
            taxa = [line.strip() for line in f if line.strip()]
        return taxa, None
    else:
        return None, re.compile(taxon_arg)

def extract_label(label: str,
                  taxon_list: Optional[List[str]],
                  taxon_regex: Optional[re.Pattern],
                  keep_accession: bool = False) -> str:
    taxon: Optional[str] = None
    if taxon_list:
        for t in taxon_list:
            if re.search(rf'(?:^|_)' + re.escape(t) + r'(?:_|$)', label):
                taxon = t
                break
    elif taxon_regex:
        m = taxon_regex.search(label)
        if m:
            taxon = m.group(0)
    if not taxon:
        return label
    if keep_accession:
        acc_m = re.match(r'^(?:[A-Z]{2,}_[0-9]+(?:\.[0-9]+)?|ENS[A-Z]*[0-9]+(?:\.[0-9]+)?)', label)
        if acc_m:
            return f"{acc_m.group(0)}_{taxon}"
    return taxon

def clean_fasta(infile: str,
                taxon_list: Optional[List[str]],
                taxon_regex: Optional[re.Pattern],
                keep_accession: bool,
                outfile: str) -> None:
    if SeqIO is None:
        sys.exit("BioPython is required for FASTA operations. Install biopython.")
    records = []
    for rec in SeqIO.parse(infile, "fasta"):
        new_id = extract_label(rec.id, taxon_list, taxon_regex, keep_accession)
        rec.id = new_id
        rec.description = ""
        records.append(rec)
    SeqIO.write(records, outfile, "fasta")
    print(f"[✓] FASTA → {outfile}")

def clean_newick(infile: str,
                 taxon_list: Optional[List[str]],
                 taxon_regex: Optional[re.Pattern],
                 keep_accession: bool,
                 outfile: str) -> None:
    if Phylo is None:
        sys.exit("BioPython is required for Newick operations. Install biopython.")
    tree = Phylo.read(infile, "newick")
    for tip in tree.get_terminals():
        tip.name = extract_label(tip.name, taxon_list, taxon_regex, keep_accession)
    Phylo.write([tree], outfile, "newick")
    print(f"[✓] Newick → {outfile}")

def run_clean_labels(args: argparse.Namespace) -> None:
    taxon_list, taxon_regex = load_taxon_filter(args.taxon_filter)
    files: List[str] = []
    for pat in args.input:
        files.extend(glob.glob(pat))
    if not files:
        sys.exit("No input files matched your patterns.")
    for fp in files:
        if not os.path.isfile(fp):
            continue
        base, ext = os.path.splitext(fp)
        out = fp if args.inplace else f"{base}_cleaned{ext}"
        try:
            if args.format == "fasta":
                clean_fasta(fp, taxon_list, taxon_regex, args.keep_accession, out)
            else:
                clean_newick(fp, taxon_list, taxon_regex, args.keep_accession, out)
        except Exception as e:
            print(f"[!] Error on {fp}: {e}", file=sys.stderr)


def register(sub, formatter_class) -> None:
    """Register the `clean-labels` subcommand."""
    sp_clean = sub.add_parser(
        "clean-labels", 
        help=("Clean FASTA or Newick labels.\n"
              "(run 'forge clean-labels -h' for options)"
        ),
    formatter_class=formatter_class,
    )    
    sp_clean.add_argument("-i", "--input", nargs="+", required=True, help="Input files (glob patterns OK)")
    sp_clean.add_argument("-f", "--format", choices=["fasta", "newick"], required=True, help="Input format")
    sp_clean.add_argument("-t", "--taxon_filter", required=True,
                          help="Regex OR path to a file listing taxa (one per line)")
    sp_clean.add_argument("-k", "--keep_accession", action="store_true", help="Keep accession prefix in labels")
    sp_clean.add_argument("--inplace", action="store_true", help="Write back to same file (no _cleaned copy)")
    sp_clean.set_defaults(func=run_clean_labels)
