module LineageCollapse
    using CSV
    using DataFrames
    using ProgressMeter
    using Clustering
    using BioSequences
    using CairoMakie

    export load_data, preprocess_data, process_lineages, plot_diagnostics

    include("data_loading.jl")
    include("preprocessing.jl")
    include("lineage_processing.jl")
    include("visualization.jl")
end
