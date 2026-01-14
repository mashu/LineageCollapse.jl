# LineageCollapse.jl
![LineageCollapse.jl Logo](https://github.com/user-attachments/assets/80776deb-d571-457e-8596-dfd0e2f834b2)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/LineageCollapse.jl/dev/)
[![Build Status](https://github.com/mashu/LineageCollapse.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/LineageCollapse.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/LineageCollapse.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/LineageCollapse.jl)

LineageCollapse.jl is a high-performance Julia package for performing soft lineage collapsing on immune repertoire sequencing data. This package provides a robust and efficient implementation of lineage collapsing, a common task in many types of immunological analyses, including antibody repertoire analysis, T cell receptor studies, and immune repertoire diversity assessment.

## Features

- Fast and memory-efficient processing of large-scale immune repertoire data
- Implements a "soft" lineage collapsing algorithm
- Supports AIRR-compliant input data formats
- Multithreaded processing for improved performance on multi-core systems

## Installation

You can install LineageCollapse.jl using Julia's package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add LineageCollapse
```

## Quick Start

## Key Concepts

**Lineage**: A cluster of sequences grouped by `v_call`, `j_call`, and `cdr3_length`,
then subdivided by the CDR3 distance clustering step (`process_lineages`). Each
cluster gets a `lineage_id`.

**Clone**: A unique combination of `d_region`, `v_call_first`, `j_call_first`, and
`cdr3` within a lineage. Clone frequency is computed per lineage.

**Collapse strategy**:
- `Hardest()` keeps exactly one clone per lineage: the clone with the highest
  `clone_frequency`.
- `Soft(cutoff)` keeps all clones within each lineage whose `clone_frequency`
  is at or above `cutoff`.

```julia
using LineageCollapse

# Load and preprocess data
df = load_data("path/to/your/data.tsv")
preprocessed_df = preprocess_data(df)

# Default is absolute mismatch distance 1 and Hardest() collapsing.
lineages = process_lineages(preprocessed_df)

# Use a mismatch threshold as a fraction of CDR3 length (e.g. 0.2 = 20%)
# This controls how sequences are clustered into lineages.
# lineages = process_lineages(preprocessed_df, 0.2)

# Or use an explicit absolute mismatch threshold (e.g. <= 1 mismatch)
# lineages = process_lineages(preprocessed_df, 1)

# Collapse identical CDR3s but only within each lineage.
# Soft(0.2) keeps clones whose (reads in clone / reads in lineage) >= 0.2.
# Hardest() keeps only the single most frequent clone per lineage.
# collapsed = collapse_lineages(lineages, Soft(0.2))
collapsed = collapse_lineages(lineages, Hardest())
```

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

## Algorithm Overview
1. **Grouping**: Sequences are grouped by `v_call`, `j_call`, and `cdr3_length`.
2. **Distance Calculation**: Pairwise Hamming distances are computed between CDR3 sequences within each group.
3. **Clustering**: Single linkage hierarchical clustering is performed on the distance matrix.
4. **Cluster Formation**: Clusters are determine by a cutoff on distances below which clusters are merged, low value means more clusters, higher value fewer clusters.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- This package was inspired by the need for standardized lineage collapsing in immune repertoire analysis.
- We thank the AIRR Community for their efforts in standardizing immune repertoire data formats.
