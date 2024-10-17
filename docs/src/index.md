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

# Assign lineages using default length Normalized Hamming distance and Hierarchical clustering
result1 = process_lineages(preprocessed_df)

# Assign lineages using length Normalized Hamming distance
# with Hierarchical clustering but different CDR3 similarity cutoff
result2 = process_lineages(preprocessed_df, 
                                    distance_metric=NormalizedHammingDistance(), 
                                    clustering_method=HierarchicalClustering(0.1))
# Collapse identical CDR3s but only within each cluster
collapsed = combine(groupby(result, [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3]), :cdr3 => length => :count)

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

