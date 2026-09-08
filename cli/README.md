# LineageCollapse CLI

`LineageCollapseCLI` wraps the Julia package as a command-line app. It reads one delimited file with `CSV.read(path, DataFrame)`, applies optional row filters, runs `preprocess_data` → `process_lineages` → `collapse_lineages`, and writes TSV (or another delimiter) to `-o` or stdout.

## Build

From the repository root:

```bash
julia --project=cli cli/build.jl --force
```

The executable is written to `build/lineage-collapse/bin/lineage-collapse`.

### GitHub Releases

When you **publish a GitHub Release** for a tag (for example after JuliaRegistrator merges and TagBot tags `v0.0.18`), the workflow [`.github/workflows/Release-cli.yml`](../.github/workflows/Release-cli.yml) builds a **Linux x86_64** tarball and uploads it to that release as `lineage-collapse-linux-x86_64-v<version>.tar.gz`. Extract it and run `lineage-collapse/bin/lineage-collapse` from inside the unpacked tree (or add `.../bin` to `PATH`). Maintainers can also run the workflow manually via **Actions → Release CLI (Linux x64) → Run workflow**; the artifact is then available under the workflow run, not attached to a release.

For broad Linux compatibility, the build script uses `cpu_target="generic"`. PackageCompiler apps still link against the libc on the build host, so build on the **oldest** glibc environment you need to support when distributing binaries.

The build enables `filter_stdlibs=true`, `include_lazy_artifacts=false`, and sysimage `sysimage_build_args=-O3` to keep the app lean and improve runtime. **Rebuild** after upgrading the CLI (`julia --project=cli cli/build.jl --force`).

### Startup time versus throughput

The first second or two after launch is Julia loading the custom sysimage; that is normal for `PackageCompiler` apps and is **not** the same as per-file processing speed. For large inputs, set **`JULIA_NUM_THREADS`** (see below). Running with **no arguments** now prints a short error and usage and exits with status 1 instead of a `MethodError`.

## Threads

Clustering uses Julia threads (`Threads.@threads` inside `LineageCollapse`). Set threads before starting the app:

```bash
JULIA_NUM_THREADS=auto build/lineage-collapse/bin/lineage-collapse input.tsv -o clonotypes.tsv
```

When developing from source:

```bash
JULIA_NUM_THREADS=auto julia --project=cli -e 'using LineageCollapseCLI; LineageCollapseCLI.main(ARGS)' -- input.tsv -o out.tsv
```

## Example data (for copy-paste)

Use any AIRR-style IgBLAST TSV with the [required columns](../README.md#input-requirements). In this repository you can point at gzipped examples under `data/` (when present), for example:

- `data/filtered-KI_100_1.tsv.gz` — filtered repertoire table  
- `data/filtered-KI_100_1-clonotypes.tsv.gz` / `data/filtered-KI_100_1-members.tsv.gz` — related tables from the same project  

Replace paths below with your own file, e.g. a full IgBLAST export:

`/home/mateusz/Mac-IgBLAST/data/ERR2856339.tsv.gz`

**Regression check:** all flags are exercised by `cli/test_options_smoke.jl` on a few-thousand-row slice of real IgBLAST output; three full-file runs (defaults, `--no-default-filters`, `--collapse soft`) were also verified on `ERR2856339.tsv.gz`.

```bash
# Fast local QA: first ~3000 rows of a huge gzipped TSV
gzip -dc /path/to/huge.tsv.gz | head -n 3001 | gzip -c > /tmp/sample.tsv.gz
julia --project=cli cli/test_options_smoke.jl /tmp/sample.tsv.gz
```

## Running without building a binary

```bash
julia --project=cli -e 'using LineageCollapseCLI; LineageCollapseCLI.main(ARGS)' -- \
  data/filtered-KI_100_1.tsv.gz -o clonotypes.tsv --quiet
```

## Defaults

With no extra flags, the tool:

- Reads tab-separated values (`--delimiter` default is `\t`).
- Applies **`V_errors<=1`** and **`J_errors<=0`** only if those columns exist (IgBLAST often uses `V_error` / `J_error`; those are accepted as aliases for the filter names `V_errors` / `J_errors`).
- Renames an existing **`count`** column to **`consensus_count`** so lineage totals are not computed from a non-additive abundance column.
- Uses **`--cdr3-differences 1`** (absolute Hamming mismatches on CDR3 within each V/J/length group).
- Uses **`--collapse hardest`** with **`--tie-breaker most-common-vdj-nt`** among maximum-frequency clones. Use `--collapse igdiscover` for IgDiscover `representative()` parity.

## Positional argument and I/O

| Option | Description |
|--------|-------------|
| `input` | Path to input TSV/CSV; gzip is supported by CSV.jl. |
| `-o`, `--output` | Output path. Use `-` for **stdout** (still tab-separated unless `--delimiter` changes it). Default is `-`. |
| `--processed-output` | Optional path for the **pre-collapse** table (same columns as after `process_lineages`, including `lineage_id`, `cluster`, frequencies). |
| `--delimiter` | Single-character field delimiter. Use the two-character token `\t` for tab (default). |

**Examples**

```bash
lineage-collapse data/filtered-KI_100_1.tsv.gz -o clonotypes.tsv
lineage-collapse data/filtered-KI_100_1.tsv.gz --processed-output lineages.tsv -o clonotypes.tsv
lineage-collapse data/filtered-KI_100_1.tsv.gz -o - | head
```

## Filtering

| Option | Description |
|--------|-------------|
| `--no-default-filters` | Do not apply the default `V_errors` / `J_errors` rules. |
| `--filter SPEC` | Extra row filter; repeatable. `SPEC` is `columnOPvalue` with `OP` in `<=`, `>=`, `==`, `!=`, `<`, `>`. Numeric columns compare as numbers; others compare as strings. |

**Examples**

```bash
# Stricter alignment (custom columns)
lineage-collapse in.tsv.gz --filter "v_identity>=0.99" --filter "j_identity>=0.99" -o out.tsv

# Productive rearrangements only (string or bool columns work via string comparison)
lineage-collapse in.tsv.gz --filter "productive==true" -o out.tsv

# Disable defaults and only use your rules
lineage-collapse in.tsv.gz --no-default-filters --filter "V_errors<=2" -o out.tsv
```

## `count` column handling

| Option | Description |
|--------|-------------|
| *(default)* | If a `count` column exists, rename it to `consensus_count` (or `--rename-count-to`). |
| `--keep-count` | Keep `count` as-is so `Hardest` / tie-breakers use summed sequence abundance. |
| `--rename-count-to NAME` | Target name for the renamed `count` column (default `consensus_count`). |

**Examples**

```bash
lineage-collapse in.tsv.gz --keep-count -o out.tsv
lineage-collapse in.tsv.gz --rename-count-to abundance -o out.tsv
```

## Preprocessing

| Option | Description |
|--------|-------------|
| `--min-d-region-length N` | Keep rows where `length(d_region) >= N` after `preprocess_data` derives `d_region`. |
| `--deduplicate` | `unique` by `sequence` (or by `sequence` + `barcode` when combined with `--use-barcode` and `barcode` exists). |
| `--use-barcode` | Meaningful only with `--deduplicate`; uses barcode when the column is present. |

**Examples**

```bash
lineage-collapse in.tsv.gz --min-d-region-length 2 -o out.tsv
lineage-collapse in.tsv.gz --deduplicate -o out.tsv
lineage-collapse in.tsv.gz --deduplicate --use-barcode -o out.tsv
```

## CDR3 clustering threshold

Exactly **one** mode applies, in this priority order:

1. `--cdr3-similarity-percent P` — minimum similarity `P` in `[0, 100]`; internally uses normalized Hamming with cutoff `1 - P/100`.
2. `--cdr3-distance-fraction F` — maximum normalized distance `F` in `[0, 1]` (normalized Hamming).
3. `--cdr3-differences K` — maximum **absolute** Hamming mismatches (non-negative integer). Default `K=1`.

`--cdr3-distance-fraction` and `--cdr3-similarity-percent` **cannot** be used together (the program errors).

**Examples**

```bash
lineage-collapse in.tsv.gz --cdr3-differences 0 -o strict.tsv
lineage-collapse in.tsv.gz --cdr3-differences 2 -o looser.tsv
lineage-collapse in.tsv.gz --cdr3-distance-fraction 0.1 -o frac.tsv
lineage-collapse in.tsv.gz --cdr3-similarity-percent 95 -o sim.tsv
```

## Linkage

| Option | Description |
|--------|-------------|
| `--linkage` | `single` (default), `complete`, or `average` — passed to hierarchical clustering in `process_lineages`. |

```bash
lineage-collapse in.tsv.gz --linkage complete -o out.tsv
```

## Collapse mode and ties (`Hardest` only)

| Option | Description |
|--------|-------------|
| `--collapse` | `hardest` (default), `soft`, or `igdiscover`. |
| `--soft-cutoff` | For `soft`, minimum clone frequency in `[0, 1]` to keep a row (default `0.3`). |
| `--tie-breaker` | For `hardest`, how to choose a representative: `most-common-vdj-nt`, `lineage-most-common-vdj-nt`, `lowest-v-shm`, `vdj-count`, `cdr3-count`, `sequence-count`, `lexicographic`, `first`, `most-naive`. |
| `--abundance` | `count` (default; use the `count` column when present, missing=1) or `rows` (each row is 1). |

#### `most-common-vdj-nt` vs `lineage-most-common-vdj-nt`

- **`most-common-vdj-nt`** (Hardest default) — among clones already at **maximum clone frequency**, pick the `vdj_nt` with the largest summed abundance.
- **`lineage-most-common-vdj-nt` / `--collapse igdiscover`** — sum abundance over **every row in the lineage**. `--collapse igdiscover` also applies IgDiscover's `n <= 2` first-row shortcut.
| `--tie-atol` | Absolute tolerance when comparing clone frequencies for `hardest` ties (default `0.0`). |

`vdj-count` and `most-naive` need a `vdj_nt` column (IgBLAST provides `VDJ_nt`, renamed automatically). `most-naive` also uses `v_identity` and `j_identity` when present.

**Examples**

```bash
lineage-collapse in.tsv.gz --collapse soft --soft-cutoff 0.25 -o soft.tsv
lineage-collapse in.tsv.gz --tie-breaker lexicographic -o out.tsv
lineage-collapse in.tsv.gz --tie-atol 0.01 -o out.tsv
```

## Logging

| Option | Description |
|--------|-------------|
| `--quiet` | Suppresses CLI progress lines **and** silences Julia logging from `LineageCollapse` during the run (`@info` from preprocessing, etc.). |

## Help

```bash
lineage-collapse --help
# or from Julia:
julia --project=cli -e 'using LineageCollapseCLI; LineageCollapseCLI.main(["--help"])'
```

## One combined example

```bash
JULIA_NUM_THREADS=auto lineage-collapse data/filtered-KI_100_1.tsv.gz \
  --filter "productive==true" \
  --cdr3-similarity-percent 92 \
  --linkage average \
  --collapse hardest \
  --tie-breaker most-common-vdj-nt \
  --processed-output lineages.tsv \
  -o clonotypes.tsv \
  --quiet
```
