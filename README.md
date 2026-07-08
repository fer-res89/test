# phylotools

Two independent command-line tools, developed together and sharing a small
internal library, for taking OrthoFinder output through to dated gene trees:

- **`forge`** — maximizes the amount of orthogroups recovered through
  OrthoFinder and curates the data for downstream analysis (QC, label
  cleaning, deduplication, BLAST-based taxon augmentation, taxid mapping).
- **`aeon`** — runs **MCMCTree** (and CodeML) at scale, with helpers for
  quartet-based gene tree bucketing, triplet tree rearrangement, Newick
  sanitization, and automated extraction of divergence-time HPDs.

They're meant to be used in sequence — `forge` produces the curated,
QC'd gene trees and FASTAs; `aeon` takes that output through batch
divergence-time estimation — but each installs and runs independently.

---

## Highlights

### forge
- **All-taxa selection**: picks orthogroups present across all species.
- **Tree QC**: flags gene trees with non-monophyletic species multiplicity (paralog check).
- **Copy matching FASTAs**: OG → FASTA matching across common extensions.
- **Post-processing options**:
  - **Deduplicate** FASTAs (keep first hit per taxon; tolerant to header formats).
  - **Clean labels** for **FASTA and/or Newick** using a taxa list (or regex). Supports keeping database accession numbers.
- **Also included**:
  - `blast-species`: append best per-species BLAST hits to input FASTAs to easily increase taxon sampling on your dataset. Supports single- or multi-database targets, parallel jobs, and resumable runs.
  - `get-taxid`: extract a seqid→taxid map from GPFF/GenBank files.

### aeon
- **Bulk MCMCTree runs**: process one FASTA, all FASTAs in a directory, or across an entire directory tree, with hand-offs to CodeML.
- **User-overridable templates**: supply your own `.ctl` template and WAG matrix file via flags, or use embedded defaults.
- **Run-dirs smart copying**: for each subdirectory, automatically copies a parent tree named `<subdir>.tre` into that subdir and runs all FASTAs there — or match per-gene trees from a `--tree-dir` by FASTA basename instead.
- **Global `--per-dir` mode**: run *any* aeon subcommand once per subdirectory (optionally in parallel, optionally targeting a specific nested directory name) without writing a wrapper script.
- **Quartet-aware gene tree bucketing**: partition gene trees by which pair of groups forms a clade, with three levels of strictness (`bipartition`, `lenient`, `strict`), and copy matching sequence/tree files into per-topology buckets.
- **Triplet rearrangements**: generate `(g1,g3)|g2` and `(g2,g3)|g1` alternatives while preserving the source tree's root split and sibling order.
- **Newick sanitization**: strip a tree to a bare topological constraint (collapse unary nodes, drop numeric internal labels, strip branch lengths), optionally merging in calibration labels from a separate calibration tree, with a batch mode that walks a whole directory tree of quartet-named subdirectories.
- **Extract 95% HPDs** from NEXUS `.tre` files, either for a fixed set of target tips or automatically per-quartet via a quartets TSV, with a batch mode for processing many quartet/data-directory combinations at once.

---

## Installation

`forge` and `aeon` live in one repository with a shared internal library
(`phylotools_common`). Each tool's third-party dependencies are gated
behind an install extra, so you only pull in what you actually need:

```bash
git clone https://github.com/<you>/phylotools.git
cd phylotools
```

> **System-wide CLI without touching system Python:** use **pipx**.

```bash
sudo apt-get update
sudo apt-get install -y pipx python3-venv
pipx ensurepath

pipx install ".[forge]"    # forge only
pipx install ".[aeon]"     # aeon only
pipx install ".[all]"      # both

# Now the CLI is available:
forge -h
aeon -h
```

> **Install inside an environment (conda/venv), e.g. for development:**

```bash
conda create -n phylotools python=3.11 -y
conda activate phylotools

python3 -m pip install -e ".[forge]"   # or .[aeon] / .[all]
```

> **Upgrading:**
```bash
git pull
pipx install ".[all]" --force            # if you used pipx
# or
python3 -m pip install -e ".[all]" --upgrade   # if you used -e .
```

**Note on the two-tool install:** `pip install .` with no extras installs
*both* console commands (`forge` and `aeon`), but neither will run until
its own extra is installed — running either one first prints a clear
message telling you which extra to add, rather than a raw traceback. This
is a structural quirk of Python packaging: a package always registers
every command it declares, regardless of which extras were selected;
extras only control which *dependencies* get pulled in.

---

## Requirements

- **Python** ≥ 3.10

- **forge** (installed via `pip install .[forge]`):
  - `biopython`, `rich` (installed automatically)
  - For `blast-species`: **BLAST+** tools (`blastn`/`blastp`/`tblastn`, `blastdbcmd`) on PATH:
    ```bash
    sudo apt update
    sudo apt install ncbi-blast+
    blastn -version
    ```
    or follow the instructions at https://www.ncbi.nlm.nih.gov/books/NBK569861
  - A BLAST database built with `-parse_seqids` for `blastdbcmd` lookups. For more details run:
    ```bash
    forge blast-species -h
    ```

- **aeon** (installed via `pip install .[aeon]`):
  - `biopython`, `dendropy`, `ete3` (installed automatically)
  - **PAML** (`mcmctree` & `codeml`) on PATH:
    - **Conda**: `conda install -c bioconda paml`
    - **Ubuntu/Debian**: `sudo apt install paml` (package availability varies by distro)
    - Or compile from source (see PAML docs)

  Verify PAML tools are on PATH:
  ```bash
  mcmctree 2>&1 | head -n 1
  codeml    2>&1 | head -n 1
  ```

---

## Commands overview

```text
forge [COMMAND] [options]
aeon  [COMMAND] [options]
```

### forge commands

| Command | Description |
|---|---|
| `filter` | Copy all-taxa trees, QC iso/paralog monophyly, copy matching FASTAs; optional dedupe and/or clean-labels. |
| `clean-labels` | Clean FASTA or Newick labels. |
| `dedupe` | Deduplicate FASTAs in-place by taxon (first match kept). |
| `blast-species` | Append best BLAST hit per species (taxid) to each FASTA. |
| `get-taxid` | Build `taxid_map.tsv` from GPFF files. |

### aeon commands

`aeon` also accepts global flags **before** the subcommand: `--per-dir`,
`--per-dir-root`, `--per-dir-include-hidden`, `--per-dir-jobs`,
`--per-dir-recursive-target` — these run the chosen subcommand once per
subdirectory (see [Per-directory execution](#per-directory-execution)
below) instead of once in the current directory.

| Command | Description |
|---|---|
| `bucket-files` | Bucket gene trees by quartet support (three modes of strictness) and copy matching sequence or tree files. |
| `process` | Run MCMCTree + CodeML on a single FASTA. |
| `pipe-all` | Run the pipeline on all FASTAs in the current directory. |
| `run-dirs` | Run the pipeline across each subdirectory (matching trees by name or via `--tree-dir`). |
| `extract-dates` | Extract 95% HPD ages from NEXUS `.tre` files, optionally per-quartet or in batch. |
| `rearrange-triplet` | Output `(g1,g3)|g2` and `(g2,g3)|g1` Newicks while preserving root split/order. |
| `sanitize-tree` | Sanitize a Newick file (collapse unary, drop numeric internal labels, strip lengths, merge calibration labels, etc.), single-file or batch. |

For full options on any command:
```bash
forge [COMMAND] -h
aeon  [COMMAND] -h
```

---

## Usage

### `forge filter`
Copy all-taxa trees, QC iso/paralog monophyly, copy matching FASTAs; optional dedupe and/or clean-labels.

**Key flags**
- `--base-dir PATH` — parent dir with the required OrthoFinder subfolders (default: `.`)
- `--trees PATH` — override the gene-trees input (a directory of per-OG tree files, or a multi-tree file). If omitted, searched for under `--base-dir`.
- `--ogseq PATH` — override the Orthogroup_Sequences directory. If omitted, searched for under `--base-dir`.
- `--gene-counts PATH` — override `Orthogroups.GeneCount.tsv`. Defaults to `--base-dir/Orthogroups/Orthogroups.GeneCount.tsv`.
- `--species-tree PATH` — override the species tree file (only required for `--dedupe` without `--taxa-file`).
- `--work-dir PATH` — where generated intermediates (e.g. split trees) are written. Default: `--base-dir/.forge_work/`
- `--dedupe` — deduplicate FASTAs in-place (first record per taxon)
- `--taxa-file FILE` — taxa list for dedupe (one per line). If omitted, falls back to `Species_Tree/SpeciesTree_rooted.txt`
- `--clean-labels TAXA_OR_REGEX` — clean labels on FASTA and Newick outputs using a taxa list file or regex
- `--keep-accession` — when cleaning, preserve accession as prefix (e.g., `XP_12345_Taxon`)
- `--inplace` — apply cleaning in-place; without it, `_cleaned` copies are created

**Outputs**
- `Gene_Trees/all_taxa_orthogroups/` — copied trees + `good_tree_files.list`, `bad_tree_files.list`
- `Orthogroup_Sequences/all_taxa_fas/` — copied FASTAs, and possibly deduped/cleaned variants
- `missing_fastas_for_good_trees.txt` (if any OGs lacked a FASTA)

---

### `forge clean-labels`
Normalize labels in FASTA or Newick using a taxa list (one per line) or a regex. Works for both trees and FASTA sequence files; supports `--keep-accession` and `--inplace`.

---

### `forge dedupe`
Deduplicate FASTAs in-place, keeping the first sequence per taxon. Species detection is tolerant: it matches normalized substrings in headers.

Can be pointed at an arbitrary FASTA directory directly with `--fasta-dir`
(default: `Orthogroup_Sequences/all_taxa_fas` under `--base-dir`), rather
than only running as part of `filter`.

---

### `forge blast-species`
For each FASTA in the current directory (matching `-e`), run BLAST and append the best hit per species (taxid) to `<name>.extended<ext>`. Produces `<name>.blast.tsv` and `<name>.extended<ext>`.

Supports single- or multi-database targets (a directory of many DBs), with
`--jobs`/`--db-threads` controlling parallelism across databases, `--extra`
for passing raw BLAST flags through, `--append {full,range,smart}` for
controlling how much of a hit gets appended, and `--resume` to skip FASTAs
whose output already indicates completion.

---

### `forge get-taxid`
Parse GPFF/GenBank files and emit a TSV mapping seqid → taxid.

---

### `aeon process`
Run MCMCTree + CodeML on one FASTA with one species tree.

```bash
aeon process OG0001.fa species_tree.tre \
  --ctl-template my_template.ctl \
  --wag-file my_wag.dat
```

- If `--ctl-template` is omitted, an embedded template is used. `seqfile=` and `treefile=` are injected automatically.
- If `--wag-file` is omitted, an embedded WAG matrix is written as `wag.dat` and used by CodeML.
- The tree can also be supplied via `--tree-dir`/`--tree-ext` instead of positionally, matched by FASTA basename.
- Pipeline: MCMCTree → patch `tmp0001.ctl` for CodeML → CodeML → write `in.BV` → final MCMCTree. Outputs are renamed alongside the FASTA (e.g., `OG0001.mcmc.txt`, `OG0001.out.txt`, `OG0001.FigTree.tre`, `OG0001.in.BV`).
- If a preliminary run fails but expected intermediates exist, the pipeline continues (with warnings).

---

### `aeon pipe-all`
Run the same analysis for every `*.fa` (or `-e` extension) in the current directory:

```bash
aeon pipe-all species_tree.tre -e .fa --resume \
  --ctl-template my_template.ctl \
  --wag-file my_wag.dat
```

- `--resume` skips FASTAs that already produced `<name>.FigTree.tre`.
- Same template, WAG, and `--tree-dir` override behavior as `process`.

---

### `aeon run-dirs`
Process each subdirectory of the current working directory.

```bash
aeon run-dirs -e .fa --resume \
  --ctl-template my_template.ctl \
  --wag-file my_wag.dat
```

**How trees are selected & copied**
- For each subdirectory `<subdir>/`, the program looks in the **parent directory** for a tree file named `<subdir>.tre`. If found, it's copied into `<subdir>/` and used for all FASTAs there.
- Alternatively, use `--tree-dir` (and `--tree-ext`) to match a tree per FASTA by basename, e.g. FASTA `OG0001.fa` + `--tree-dir trees/` → `trees/OG0001.tre`.
- With `--tree-dir-search`, when processing subdirectory `X`, trees are searched for under `<tree-dir>/X/` instead, matched by FASTA basename; if that directory doesn't exist, `X` is skipped.
- `--sanitize` runs `sanitize-tree` on per-gene trees pulled from `--tree-dir` before use, using the same letter-option syntax as `sanitize-tree` (see below).

**Example layout**
```
parent_dir/
  ├── pal_neog/
  │   ├── OG0001.fa
  │   └── OG0002.fa
  ├── pal_gall/
  │   └── OG0003.fa
  ├── pal_neog.tre    # used for FASTAs in pal_neog/
  └── pal_gall.tre    # used for FASTAs in pal_gall/
```

---

### `aeon extract-dates`
Extract 95% HPD bounds (and their average) for a clade defined by target tips:

```bash
aeon extract-dates --target-tips TaxonA,TaxonB -e .tre
# writes: dates_<cwd>.csv (or dates_int_/dates_ext_ with -in/-ex)
```

Parses NEXUS trees, finds the common ancestor of the provided tips, and reads the `&95%HPD={lo,hi}` comment. The CSV has columns: `base, lo, hi, avg`.

Instead of fixed `--target-tips`, you can drive this from a quartets TSV
(same format as `bucket-files`) with `-qt`/`--quartets-tsv` plus
`-auto`/`--auto-quartets`: for each quartet directory found, `A_set∪B_set`
is used as the internal-date ingroup and `A_set∪B_set∪C_set` as the
external-date set. `-out`/`--out-dir` controls where per-quartet,
per-data-directory CSVs land. `-b`/`--batch` extends this further: starting
from the current directory, it looks for subdirectories matching a
quartet's `name` and processes each of *their* immediate subdirectories as
data directories.

---

### `aeon bucket-files`
Bucket gene trees by which pair of groups forms a clade (G1+G2, G1+G3, G2+G3), and copy matching sequence or tree files.

```bash
# groups.txt has exactly 3 lines: group1, group2, group3 (comma-separated taxa)
aeon bucket-files -e .tre -G groups.txt -f Orthogroup_Sequences/all_taxa_fas -s .fa \
  -o G1G2_bucket,G1G3_bucket,G2G3_bucket
```

Creates output directories and copies the corresponding files whose
basenames match the input tree files. `--bucket-mode` controls how
strictly quartet support is interpreted:
- `bipartition` (default) — unrooted, split-based (MSC-style) classification.
- `lenient` — bipartition classification, then requires some rooted constraints (focal clades roughly monophyletic) for an informative label.
- `strict` — as `lenient`, with stronger rooted-clade requirements; `--max-strays` caps how many focal-taxon strays are tolerated before a tree is demoted to a `root_fail` bucket.

`--bs-thresh` collapses low-support edges before classification;
`--min-per-side` sets the minimum taxa required per quartet side.

---

### `aeon rearrange-triplet`
Produce two rearranged trees from one input, keeping the original root split and sibling order. Singleton parentheses are collapsed for readability.

```bash
# taxa.txt has 3 non-empty lines: group1, group2, group3 (comma-separated)
aeon rearrange-triplet -i input.tre -t taxa.txt \
  -o g1g3_sis_g2 g2g3_sis_g1
# writes: g1g3_sis_g2.tre and g2g3_sis_g1.tre
```

`--sanitize` can sanitize the input and/or output trees inline using the
same letter-option syntax as `sanitize-tree`.

---

### `aeon sanitize-tree`
Strips a Newick file to bare bones to make it usable as a topological constraint (collapse unary nodes, drop numeric internal labels, strip lengths, etc.).

```bash
aeon sanitize-tree -i input.tre -o output.txt -s -p -H
# keeps support values (-s), allows polytomies (-p), adds a Phylip-style header (-H)
```

- By default, polytomies stop the sanitize process with a warning; allow them with `-p`/`--allow-polytomies`.
- `-H`/`--add-header` adds a `<num_taxa> <num_trees>` header, required by MCMCTree.
- `-l`/`--keep-lengths`, `-s`/`--keep-support`, `-c`/`--no-collapse` control what's stripped.
- `--cal-tree PATH` merges calibration labels from a separate calibration tree in **after** sanitizing; `--label-prefix` can prefix merged labels (e.g. `cal:`); `--no-overwrite` avoids clobbering existing non-empty internal labels.
- `--quartets-tsv`/`--quartet` supply A/B/C/D group definitions (same format as `bucket-files`) for MRCA-based label scrubbing.
- `--batch` treats `-i`/`-o` as an extension/suffix pair and walks subdirectories whose names match quartets in `--quartets-tsv`, rooted at `--root` (default: cwd); `--tree-ext` and `--suffix` override the extension/suffix independently if needed. `--dry-run` shows what would run without writing anything.

---

## Per-directory execution

Any `aeon` subcommand can be run once per subdirectory using global flags
placed **before** the subcommand:

```bash
aeon --per-dir --per-dir-root . --per-dir-jobs 4 pipe-all species_tree.tre -e .fa
```

- `--per-dir-root DIR` — root whose immediate subdirectories are processed (default: `.`).
- `--per-dir-include-hidden` — include dot-directories.
- `--per-dir-jobs N` — run up to `N` subdirectories in parallel (default: 1, sequential). Failures in one subdirectory are reported without stopping the others.
- `--per-dir-recursive-target NAME` — instead of the immediate subdirectories themselves, search inside each of them for subdirectories named exactly `NAME` (at any depth) and use those as the targets. E.g. with `--per-dir-root .` and `--per-dir-recursive-target runs`, the command runs inside every directory named `runs` found under each immediate child of `.`.

Inside each run, relative paths in other options are interpreted relative
to that subdirectory.

---

## Tips & Notes

### forge
- **Run location**: `forge filter` should be run from the parent directory of OrthoFinder results containing `Orthogroups/`, `Gene_Trees/`, `Orthogroup_Sequences/`, and `Species_Tree/` — or a directory structure that resembles that, via `--trees`/`--ogseq`/`--species-tree`/`--gene-counts`.
- **All-taxa selection**: Orthogroups are selected from `Orthogroups/Orthogroups.GeneCount.tsv` where every species is present (no zeros/NA).
- **Tree QC**: Good/bad lists are written as `Gene_Trees/all_taxa_orthogroups/good_tree_files.list` and `bad_tree_files.list`. FASTA copy uses only the good OG IDs.
- **Missing FASTAs**: If a good tree has no matching FASTA, see `Orthogroup_Sequences/all_taxa_fas/missing_fastas_for_good_trees.txt`. Check filename stems (e.g., `OG0001234.fasta`) and extensions (`.fa|.fas|.fasta`).
- **Dedupe logic**: Keeps the first sequence per taxon. Matching is tolerant (substring on normalized headers), but headers must contain the taxon label somewhere. Use `--taxa-file` to control species names explicitly; otherwise falls back to `Species_Tree/SpeciesTree_rooted.txt`.
- **Clean labels**: accepts either a taxa list file (one per line) or a regex. Works for FASTA and Newick. Add `--keep-accession` to preserve an accession prefix (e.g., `XP_12345_Taxon`).
- **In-place vs copies**: Add `--inplace` to overwrite outputs. Without it, `_cleaned` files are created alongside originals.
- **Suffix hygiene**: Inputs with `_clean` or `_edit` in the stem are skipped by dedupe to avoid reprocessing artifacts.
- **Line endings** (WSL/Windows): if editing under `/mnt/c`, convert files before running:
  ```bash
  dos2unix -k -o $(git ls-files '*.py' '*.fa*' '*.txt' 2>/dev/null)
  ```

### aeon
- **Templates**: when providing a custom `.ctl` via `--ctl-template`, `seqfile=`/`treefile=` lines are patched to the current FASTA/tree names automatically.
- **Matrix files**: `--wag-file` copies your file to `wag.dat` in the working directory so CodeML finds it via `aaRatefile = wag.dat`.
- **Taxon labels**: underscores are preserved and taxa are matched literally — ensure FASTA and tree labels agree.
- **HPD extraction**: requires NEXUS trees containing `&95%HPD={...}` annotations on clade comments.
- **Resume mode**: `--resume` is conservative; it skips entries if `<name>.FigTree.tre` exists. Delete or move that file to force a re-run.
- **Header on tree file**: `mcmctree` only accepts trees that include a header with the number of taxa and trees included in the file — use `sanitize-tree -H` to generate one.
- **Unary nodes**: `rearrange-triplet` can sometimes produce unary nodes; sanitize the resulting trees afterward if so.
- **Quartet-driven workflows**: `bucket-files`, `extract-dates -auto`, and `sanitize-tree --quartets-tsv`/`--batch` all consume the same quartets TSV format (`name, A_set, B_set, C_set, D_set`), so one quartet definition file can drive bucketing, tree sanitization, and date extraction across the whole pipeline.

---

## Troubleshooting

- **`forge`/`aeon`: command not found**
  Run `pipx ensurepath`, open a new terminal, or `source ~/.bashrc`. Verify install with `pipx list`.

- **Debian "externally managed environment" error**
  Use `pipx install --editable ".[forge]"` (or `.[aeon]`/`.[all]`), or a user install: `python3 -m pip install --user -e ".[forge]"`. Avoid system-wide pip unless you use `--break-system-packages`.

- **`forge is missing a required dependency (...)` / `aeon is missing a required dependency (...)`**
  You installed the package without the matching extra. Run `pip install "phylotools[forge]"` or `pip install "phylotools[aeon]"` (or `[all]`) in the same environment. If using pipx: `pipx runpip forge install biopython rich` / `pipx runpip aeon install dendropy ete3`.

- **Fewer FASTAs than good trees (forge)**
  Check `Orthogroup_Sequences/all_taxa_fas/missing_fastas_for_good_trees.txt`. Confirm filenames (e.g., `OG0001234.fas`) exactly match OG stems; verify extensions are `.fa|.fas|.fasta` (case-insensitive).

- **Dedupe "does nothing" (forge)**
  Ensure headers contain the taxon string you expect (post-cleaning they should). Try running `clean-labels` first (with your taxa list), then `dedupe`. You can also pass `--taxa-file` to align species names precisely.

- **Label cleaning "no effect" (forge)**
  Verify `--clean-labels` points to an existing taxa file (or provide a regex). Add `--keep-accession` if you want to retain accession prefixes. Without `--inplace`, look for `_cleaned` outputs.

- **Newick parse errors / empty species tree (forge)**
  Confirm `Species_Tree/SpeciesTree_rooted.txt` is valid Newick and that terminal names are present. If not available, supply a taxa list via `--taxa-file` for dedupe/cleaning.

- **`mcmctree`/`codeml` not found (aeon)**
  Make sure PAML is installed and on PATH. With conda: `conda install -c bioconda paml`.

- **Label mismatches / missing taxa (aeon)**
  Confirm FASTA headers and tree tip labels match exactly (including underscores).

- **`mcmctree` rejects a tree file (aeon)**
  It needs a header with the taxon/tree count — run it through `aeon sanitize-tree -H` first.


---

## Citation

If your work uses **forge** and/or **aeon**, please cite the tool(s) and their core dependencies.

**forge / aeon**

'ILS paper citation'

**Also cite, as appropriate:**
- **OrthoFinder** (for orthogroup inference and outputs curated with `forge`)
- **Biopython** (FASTA/Newick parsing and writing)
- **BLAST+** (if you use `forge blast-species`)
- **PAML / MCMCTree** (for `aeon`'s divergence-time estimation)
- **DendroPy** and **ETE Toolkit** (tree manipulation within `aeon`)
