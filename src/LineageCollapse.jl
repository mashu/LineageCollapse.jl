module LineageCollapse
    using CSV
    using DataFrames
    using Clustering
    using BioSequences
    using StringDistances

    export load_data, preprocess_data, deduplicate_data, process_lineages, collapse_lineages
    export hardest_tie_summary
    export AbstractDistanceMetric, ClusteringMethod, CollapseStrategy, AbstractTieBreaker, TieBreaker
    export HammingDistance, NormalizedHammingDistance, LevenshteinDistance, HierarchicalClustering
    export Hardest, Soft
    export ByVdjCount, ByCdr3Count, BySequenceCount, ByLexicographic, ByFirst, ByMostNaive, ByMostCommonVdjNt
    export compute_distance, compute_pairwise_distance, perform_clustering

    include("data_loading.jl")
    include("preprocessing.jl")
    include("tie_breaking.jl")
    include("lineage_processing.jl")
end
