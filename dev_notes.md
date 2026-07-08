## Developer notes

### Repository layout

```
src/
  phylotools_common/   # generic, dependency-light file/path helpers
                        # shared by both tools (NOT tree-parsing logic —
                        # forge uses Bio.Phylo, aeon uses dendropy/ete3,
                        # and those are kept separate)
  forge/
    cli.py              # builds the argparse parser, registers subcommands
    _common.py          # forge-specific constants + Bio.Phylo helpers
    filter_pipeline.py  # `forge filter`
    clean_labels.py     # `forge clean-labels`
    dedupe.py           # `forge dedupe`
    blast_species.py    # `forge blast-species`
    taxid.py            # `forge get-taxid`
  aeon/
    cli.py              # builds the argparse parser, registers subcommands,
                         # implements the global --per-dir wrapper
    _common.py           # embedded CTL/WAG templates, mcmctree/codeml helpers
    bucketing.py         # `aeon bucket-files`
    process.py           # `aeon process` / `pipe-all` / `run-dirs`
    extract_dates.py     # `aeon extract-dates`
    rearrange_triplet.py # `aeon rearrange-triplet`
    sanitize_tree.py     # `aeon sanitize-tree` (also imported by
                          # process.py and rearrange_triplet.py)
```

### Adding a new subcommand

Each subcommand module exposes a `register(sub[, formatter_class])`
function that adds its own argparse subparser. To add a new one:

1. Create a new module in `forge/` or `aeon/` with a `register` function
   following the pattern in any existing module, plus whatever
   `run_*`/`cmd_*` function it points to via `set_defaults(func=...)`.
2. Add the module to the `SUBCOMMAND_MODULES` tuple at the top of the
   corresponding `cli.py`.

No other file needs to change.

### Known follow-ups worth deciding on before publication

- No automated tests yet (`tests/` is a placeholder)
---

