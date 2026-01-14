```@meta
CurrentModule = LineageCollapse
```

# LineageCollapse

Documentation for [LineageCollapse](https://github.com/mashu/LineageCollapse.jl), a Julia package for collapsing lineages in AIRR data.

## Overview

LineageCollapse provides tools for processing and analyzing adaptive immune receptor repertoire (AIRR) data. It offers functions for data loading, preprocessing, lineage collapsing.

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

# Assign lineages using default absolute mismatch distance 1
result1 = process_lineages(preprocessed_df)

# Assign lineages using a CDR3 mismatch fraction (0.1 = 10%)
result2 = process_lineages(preprocessed_df, 0.1)

# Assign lineages using an absolute mismatch threshold (<= 1 mismatch)
result3 = process_lineages(preprocessed_df, 1)
# Collapse
collapsed_df = collapse_lineages(lineages, Soft(0.2))

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

