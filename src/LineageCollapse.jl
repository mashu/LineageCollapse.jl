module LineageCollapse
    using CSV
    using DataFrames
    using ProgressMeter
    using Clustering
    using BioSequences
    using StringDistances

    export load_data, preprocess_data, process_lineages, plot_diagnostics
    export DistanceMetric, ClusteringMethod
    export HammingDistance, LevenshteinDistance, HierarchicalClustering
    export compute_distance, compute_pairwise_distance, perform_clustering

    function plot_diagnostics(args...; opts...)
        if !isdefined(@__MODULE__, :CairoMakie)
            error("Import CairoMakie to enable plotting")
        else
            error("Invalid method call")
        end
    end
    include("data_loading.jl")
    include("preprocessing.jl")
    include("lineage_processing.jl")
end
