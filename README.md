# LineageCollapse.jl
![LineageCollapse.jl Logo](https://github.com/user-attachments/assets/80776deb-d571-457e-8596-dfd0e2f834b2)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/LineageCollapse.jl/dev/)
[![Build Status](https://github.com/mashu/LineageCollapse.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/LineageCollapse.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/LineageCollapse.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/LineageCollapse.jl)

LineageCollapse.jl is a high-performance Julia package for performing soft lineage collapsing on immune repertoire sequencing data. This package provides a robust and efficient implementation of lineage collapsing, a common task in many types of immunological analyses, including antibody repertoire analysis, T cell receptor studies, and immune repertoire diversity assessment.

## Features

- Fast and memory-efficient processing of large-scale immune repertoire data
- Multiple collapse strategies: `Hardest()` (one representative per lineage) and `Soft(cutoff)` (frequency threshold)
- Configurable distance metrics: Hamming, Normalized Hamming, Levenshtein
- Extensible tie-breaking system for representative selection
- Supports AIRR-compliant input data formats
- Multithreaded processing for improved performance on multi-core systems

## Installation

You can install LineageCollapse.jl using Julia's package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add LineageCollapse
```

## Quick Start

```julia
using LineageCollapse

# Load and preprocess data
df = load_data("path/to/your/data.tsv")
preprocessed_df = preprocess_data(df; min_d_region_length=3)

# Assign lineages using absolute mismatch threshold (≤ 1 mismatch)
lineages = process_lineages(preprocessed_df, 1)

# Or use a mismatch fraction of CDR3 length (0.1 = 10%)
# lineages = process_lineages(preprocessed_df, 0.1)

# Collapse to one representative per lineage (Hardest)
collapsed = collapse_lineages(lineages, Hardest())

# Or keep clones with frequency ≥ 20% (Soft)
# collapsed = collapse_lineages(lineages, Soft(0.2))
```

## Key Concepts

| Term | Description |
|------|-------------|
| **Lineage** | Sequences grouped by V-gene, J-gene, CDR3 length, and clustered by CDR3 similarity |
| **Clone** | Unique combination of D-region + V + J + CDR3 within a lineage |
| **Clone Frequency** | Proportion of sequences in a lineage belonging to a clone |
| **Hardest** | Collapse strategy: one representative per lineage (highest frequency) |
| **Soft** | Collapse strategy: keep all clones above a frequency threshold |

## Input Requirements

LineageCollapse.jl requires input data to be in AIRR-C (Adaptive Immune Receptor Repertoire - Community) format, typically obtained from tools like IgBLAST. The following columns are required:

- sequence_id
- sequence
- v_sequence_end
- j_sequence_start
- cdr3
- v_call
- j_call
- stop_codon

## Tie-Breaking

When multiple clones have the same maximum frequency in `Hardest()` mode, a tie-breaker selects the representative:

```julia
# Default: most common VDJ nucleotide sequence (igdiscover-compatible)
collapsed = collapse_lineages(lineages, Hardest(); tie_breaker=ByMostCommonVdjNt())

# Alternative: prioritize sequences closest to germline
collapsed = collapse_lineages(lineages, Hardest(); tie_breaker=ByMostNaive())

# With frequency tolerance (1% tolerance for "tied" clones)
collapsed = collapse_lineages(lineages, Hardest(); tie_atol=0.01)
```

Available tie-breakers: `ByMostCommonVdjNt()`, `ByVdjCount()`, `ByCdr3Count()`, `BySequenceCount()`, `ByMostNaive()`, `ByLexicographic()`, `ByFirst()`

## Algorithm Overview

1. **Grouping**: Sequences are grouped by `v_call`, `j_call`, and `cdr3_length`.
2. **Distance Calculation**: Pairwise distances are computed between CDR3 sequences within each group (Hamming, Normalized Hamming, or Levenshtein).
3. **Clustering**: Hierarchical clustering (default: single linkage) is performed on the distance matrix.
4. **Cluster Formation**: Clusters are formed by cutting the dendrogram at the specified threshold. Lower values create more clusters, higher values fewer clusters.
5. **Collapsing**: Representatives are selected per lineage (`Hardest`) or clones above a frequency threshold are retained (`Soft`).

## Benchmark

Performance on a 617K sequence repertoire (Intel Core i9-13980HX, 32 threads):

| Stage | Time | % of Total | Rate |
|-------|------|------------|------|
| IO (gzipped TSV) | 6.2 s | 64.6% | 100K seq/s |
| Preprocessing | 1.2 s | 13.1% | 496K seq/s |
| Clustering | 0.4 s | 4.3% | 1.5M seq/s |
| Collapse (Hardest) | 1.7 s | 18.1% | 358K seq/s |
| **Total** | **9.5 s** | 100% | **65K seq/s** |

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- This package was inspired by the need for standardized lineage collapsing in immune repertoire analysis.
- We thank the AIRR Community for their efforts in standardizing immune repertoire data formats.
