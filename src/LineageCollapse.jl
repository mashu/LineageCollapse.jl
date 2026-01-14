module LineageCollapse
    using CSV
    using DataFrames
    using ProgressMeter
    using Clustering
    using BioSequences
    using StringDistances
    using SparseArrays
    using LinearAlgebra

    export load_data, preprocess_data, deduplicate_data, process_lineages, collapse_lineages
    export DistanceMetric, ClusteringMethod, CollapseStrategy
    export HammingDistance, NormalizedHammingDistance, LevenshteinDistance, HierarchicalClustering
    export Hardest, Soft
    export compute_distance, compute_pairwise_distance, perform_clustering

    include("data_loading.jl")
    include("preprocessing.jl")
    include("lineage_processing.jl")
end
