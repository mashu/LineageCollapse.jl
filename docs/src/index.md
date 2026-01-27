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
| `ByMostCommonVdjNt()` | Most common VDJ nucleotide sequence (igdiscover-compatible) |
| `ByVdjCount()` | Highest VDJ count, then lexicographic CDR3 |
| `ByCdr3Count()` | Highest CDR3 count, then lexicographic CDR3 |
| `BySequenceCount()` | Highest sequence count, then lexicographic CDR3 |
| `ByMostNaive()` | Highest V/J identity (closest to germline) |
| `ByLexicographic()` | Lexicographically smallest CDR3 |
| `ByFirst()` | First candidate (no sorting) |

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

For detailed function signatures and options, see the API Reference below.

```@index
```

```@autodocs
Modules = [LineageCollapse]
```

