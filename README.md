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

```julia
using LineageCollapse

# Load and preprocess data
df = load_data("path/to/your/data.tsv")
preprocessed_df = preprocess_data(df)

# Use default Normalized Hamming distance with Hierarchical clustering and CDR3 similarity cutoff of 0.2
result = process_lineages(preprocessed_df, clustering_method=HierarchicalClustering(0.2))

# Collapse identical CDR3s but only within each cluster
# collapsed = collapse_lineages(lineages, 0.2, :soft)
collapsed = collapse_lineages(lineages, 0.2, :hardrest)
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
