```@meta
CurrentModule = LineageCollapse
```

# LineageCollapse

Documentation for [LineageCollapse](https://github.com/mashu/LineageCollapse.jl), a Julia package for collapsing lineages in AIRR data.

## Overview

LineageCollapse provides tools for processing and analyzing adaptive immune receptor repertoire (AIRR) data. It offers functions for data loading, preprocessing, lineage assignment, and lineage collapsing.

## Architecture

### Processing Pipeline

The library follows a linear pipeline where each stage transforms your data:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                          LineageCollapse Pipeline                           │
└─────────────────────────────────────────────────────────────────────────────┘

   AIRR TSV File
        │
        ▼
┌───────────────┐
│  load_data()  │  Load and validate AIRR-formatted data
└───────┬───────┘
        │ DataFrame
        ▼
┌──────────────────────┐
│  preprocess_data()   │  Filter, derive D-region, normalize columns
└──────────┬───────────┘
           │ DataFrame + d_region, v_call_first, j_call_first
           ▼
┌────────────────────────────────────────────────────────────────────────────┐
│  process_lineages(df, threshold)                                           │
│  ┌──────────────────────────────────────────────────────────────────────┐  │
│  │  For each (V-gene, J-gene, CDR3-length) group:                       │  │
│  │    1. Compute pairwise CDR3 distances (AbstractDistanceMetric)       │  │
│  │    2. Cluster sequences (ClusteringMethod)                           │  │
│  │    3. Assign lineage_id to each cluster                              │  │
│  └──────────────────────────────────────────────────────────────────────┘  │
└──────────┬─────────────────────────────────────────────────────────────────┘
           │ DataFrame + lineage_id, cluster, cluster_size, cdr3_frequency
           ▼
┌────────────────────────────────────────────────────────────────────────────┐
│  collapse_lineages(df, strategy)                                           │
│  ┌──────────────────────────────────────────────────────────────────────┐  │
│  │  CollapseStrategy determines output:                                 │  │
│  │    • Hardest() → one representative per lineage                      │  │
│  │    • Soft(cutoff) → all clones above frequency threshold             │  │
│  │                                                                      │  │
│  │  AbstractTieBreaker resolves ties when selecting representatives     │  │
│  └──────────────────────────────────────────────────────────────────────┘  │
└──────────┬─────────────────────────────────────────────────────────────────┘
           │
           ▼
   Collapsed DataFrame (representatives or filtered clones)
```

### Type Hierarchy

The library uses Julia's multiple dispatch with abstract types for extensibility:

```
AbstractDistanceMetric           ClusteringMethod              CollapseStrategy
         │                              │                             │
         ├── HammingDistance            └── HierarchicalClustering    ├── Hardest
         ├── NormalizedHammingDistance          │                     └── Soft
         └── LevenshteinDistance                └── cutoff::Float
                                                                           │
                                                                           │
AbstractTieBreaker ◄───────────────────────────────────────────────────────┘
         │                                                    (used by Hardest)
         ├── TieBreaker
         │       └── criteria::Vector{Pair{Symbol,Bool}}
         │
         └── MostCommonVdjNtTieBreaker
```

### Extension Points

You can extend the library by defining new types and methods:

**Custom Distance Metric:**
```julia
struct MyDistance <: AbstractDistanceMetric end

# Implement the required method
LineageCollapse.compute_distance(::MyDistance, x::LongDNA{4}, y::LongDNA{4}) = ...
```

**Custom Tie-Breaker:**
```julia
struct MyTieBreaker <: AbstractTieBreaker end

# Or use the built-in TieBreaker with custom criteria
my_breaker = TieBreaker([:my_column => true, :cdr3 => false])
```

### Key Concepts

| Concept | Description |
|---------|-------------|
| **Lineage** | Group of sequences sharing V-gene, J-gene, CDR3 length, and similar CDR3 |
| **Clone** | Unique combination of D-region + lineage + V + J + CDR3 within a lineage |
| **Clone Frequency** | Proportion of sequences in a lineage belonging to a clone |
| **Representative** | Selected sequence to represent an entire lineage |

### Built-in Tie-Breakers

| Function | Strategy |
|----------|----------|
| `ByMostCommonVdjNt()` | Most common `vdj_nt` by **lineage-wide** summed `count` (igdiscover-style); see note below vs `ByVdjCount` |
| `ByVdjCount()` | Among **max-frequency** candidates, highest **per-clone** `vdj_count`, then lexicographic `cdr3` |
| `ByCdr3Count()` | Highest CDR3 count, then lexicographic CDR3 |
| `BySequenceCount()` | Highest sequence count, then lexicographic CDR3 |
| `ByMostNaive()` | Highest V/J identity (closest to germline) |
| `ByLexicographic()` | Lexicographically smallest CDR3 |
| `ByFirst()` | First candidate (no sorting) |

**`ByMostCommonVdjNt()` vs `ByVdjCount()`:** `ByMostCommonVdjNt()` pools every row in the lineage to decide which `vdj_nt` has the largest total abundance, then picks the first row with that sequence (that row need not be on a maximum-frequency clone). `ByVdjCount()` only compares rows already tied at **maximum clone frequency** and ranks them using a score computed **inside each clone**: the maximum, over distinct `vdj_nt` in the clone, of the summed `count` for that `vdj_nt`. Use the first for igdiscover parity; use the second when the representative must be drawn from the top-frequency clone set with a local VDJ-based score.

Tie-breakers can be combined: `ByVdjCount() + ByLexicographic()`

## Installation

```julia
using Pkg
Pkg.add("LineageCollapse")
```

## Quick Start

```julia
using LineageCollapse

# Load and preprocess
df = load_data("airr_data.tsv.gz")
df = preprocess_data(df; min_d_region_length=3)

# Assign lineages (threshold: 1 mismatch or 0.1 = 10% of CDR3 length)
lineages = process_lineages(df, 1)

# Collapse options:
# Option A: One representative per lineage
result = collapse_lineages(lineages, Hardest())

# Option B: Keep clones with frequency ≥ 20%
result = collapse_lineages(lineages, Soft(0.2))

# Option C: Custom tie-breaking
result = collapse_lineages(lineages, Hardest(); 
                           tie_breaker=ByMostNaive(),
                           tie_atol=0.01)  # 1% tolerance
```

## Detailed Examples

### Using Different Distance Metrics

```julia
# Explicit metric configuration
lineages = process_lineages(df;
    distance_metric = NormalizedHammingDistance(),
    clustering_method = HierarchicalClustering(0.1f0),
    linkage = :average
)
```

### Diagnostic: Finding Ties

```julia
# Identify lineages where multiple clones have the same max frequency
ties = hardest_tie_summary(df; atol=0.01)
filter(:hardest_tied => identity, ties)  # Show only tied lineages
```

## Column reference

Meaning and computation of columns added by each stage. The most important is **`count`** (optional input and main abundance output).

### Count (most important)

**Input (optional):** An optional numeric column `count` (e.g. read or UMI count per sequence). If present, it is the *abundance* of that row; if absent, each row is treated as abundance 1. Missing values are treated as 0 when summing.

**Output (Hardest only):** The collapsed result includes a `count` column: for each lineage it is the **sum of the input `count`** over all sequences in that lineage (total lineage abundance). If the input had no `count` column, this equals the number of sequences in the lineage.

**How it is used:** `ByMostCommonVdjNt()` and `ByVdjCount()` use summed `count` when ranking `vdj_nt` choices (see **Built-in Tie-Breakers** above for how they differ). Clone frequency and the Soft strategy use the **number of sequences (rows)** per clone, not the sum of `count`.

### After `preprocess_data`

| Column | Meaning | Computation |
|--------|---------|-------------|
| `d_region` | D-region nucleotide sequence | `sequence[v_sequence_end+1:j_sequence_start]` |
| `v_call_first` | First V-gene allele | First token of `v_call` (before comma) |
| `j_call_first` | First J-gene allele | First token of `j_call` (before comma) |
| `vdj_nt` | V–D–J nucleotide sequence | `sequence[v_sequence_start:j_sequence_end]` (only if those columns exist) |
| `cdr3_length` | CDR3 length | `length(cdr3)` |

### After `process_lineages`

| Column | Meaning | Computation |
|--------|---------|-------------|
| `lineage_id` | Lineage identifier | Unique id per (V, J, CDR3 length, cluster) |
| `cluster` | Cluster within V/J/CDR3-length group | From hierarchical clustering of CDR3 distances |
| `cluster_size` | Sequences in this cluster | Number of rows in the same cluster |
| `min_distance` | Min distance to another CDR3 in cluster | Smallest pairwise distance to another sequence in the cluster |
| `cdr3_count` | How many sequences share this CDR3 in the cluster | Number of rows with same CDR3 in same cluster |
| `max_cdr3_count` | Max `cdr3_count` in this cluster | Maximum of `cdr3_count` over the cluster |
| `cdr3_frequency` | Relative frequency of this CDR3 in cluster | `cdr3_count / max_cdr3_count` (0–1) |

### After `collapse_lineages` with Hardest

| Column | Meaning | Computation |
|--------|---------|-------------|
| `count` | **Total lineage abundance** | Sum of input `count` (or 1 per row if no input `count`) over all sequences in the lineage — see "Count (most important)" above |
| `nVDJ_nt` | Number of unique VDJ nucleotide sequences in lineage | Count of distinct `vdj_nt` in the lineage (or `missing` if no `vdj_nt`) |

### After `collapse_lineages` with Soft

| Column | Meaning | Computation |
|--------|---------|-------------|
| `clone_frequency` | Fraction of lineage (by row count) in this clone | Number of sequences (rows) in this clone ÷ total sequences in the lineage |
| `sequence_count` | Number of sequences in this clone | Number of rows with same (d_region, lineage_id, v_call_first, j_call_first, cdr3) |

## Command-line interface

A standalone CLI (`lineage-collapse`) lives under [`cli/`](https://github.com/mashu/LineageCollapse.jl/tree/main/cli) in the source tree. It reads one TSV (or gzip-compressed TSV), applies optional column filters, runs the same preprocessing and lineage pipeline as the Julia API, and writes collapsed output.

**Build a binary** (from a clone of the repository):

```bash
julia --project=cli cli/build.jl --force
```

The app is created under `build/lineage-collapse/`; the executable is `build/lineage-collapse/bin/lineage-collapse`.

**Threads:** set `JULIA_NUM_THREADS` (e.g. `auto`) before launching the binary so `process_lineages` can use multiple threads.

**Run without building** (same arguments as the binary):

```bash
julia --project=cli -e 'using LineageCollapseCLI; LineageCollapseCLI.main(ARGS)' -- path/to/input.tsv.gz -o out.tsv
```

**Reference and examples:** every flag, defaults, example data paths, and usage patterns are documented in [cli/README.md](https://github.com/mashu/LineageCollapse.jl/blob/main/cli/README.md) in the repository.

**Linux binary:** each published GitHub Release can include a prebuilt `lineage-collapse` app tarball built by [`.github/workflows/Release-cli.yml`](https://github.com/mashu/LineageCollapse.jl/blob/main/.github/workflows/Release-cli.yml) (Linux x86_64). Download the asset from the release page, unpack, and run `bin/lineage-collapse`.

**Regression driver** for CLI options (pass any AIRR-style TSV or `.gz`):

```bash
julia --project=cli cli/test_options_smoke.jl path/to/example.tsv.gz
```

For detailed function signatures and options of the Julia API, see the API Reference below.

```@index
```

```@autodocs
Modules = [LineageCollapse]
```

