```@meta
CurrentModule = LineageCollapse
```

# LineageCollapse

Documentation for [LineageCollapse](https://github.com/mashu/LineageCollapse.jl), a Julia package for collapsing lineages in AIRR data.

## Overview

LineageCollapse provides tools for processing and analyzing adaptive immune receptor repertoire (AIRR) data. It offers functions for data loading, preprocessing, lineage collapsing, and visualization.

## Installation

You can install LineageCollapse using Julia's package manager:

```julia
using Pkg
Pkg.add("LineageCollapse")
```

## Basic Usage

Here's a quick example of how to use LineageCollapse:

```julia
using LineageCollapse

# Load data
df = load_data("path/to/your/airr_data.tsv.gz")

# Preprocess data
preprocessed_df = preprocess_data(df, min_d_region_length=3)

# Perform lineage collapsing using default Hamming distance and Hierarchical clustering
collapsed_df = process_lineages(preprocessed_df)

# Use Levenshtein distance with Hierarchical clustering
collapsed_df_lev = process_lineages(preprocessed_df, 
                                    distance_metric=LevenshteinDistance(), 
                                    clustering_method=HierarchicalClustering(0.1))

# Adjust allele ratio and collapse results
collapsed_df_custom = process_lineages(preprocessed_df, 
                                       distance_metric=NormalizedHammingDistance(),
                                       clustering_method=HierarchicalClustering(0.1),
                                       cdr3_ratio=0.3,
                                       collapse=true)

# Generate diagnostic plots (requires CairoMakie)
# using CairoMakie
# plot_diagnostics(collapsed_df)
```

For more detailed information on each function and its options, please refer to the API documentation below.

```@index
```

```@autodocs
Modules = [LineageCollapse]
```

