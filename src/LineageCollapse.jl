module LineageCollapse
    using CSV
    using DataFrames
    using ProgressMeter
    using Clustering
    using BioSequences

    export load_data, preprocess_data, process_lineages, plot_diagnostics
    function plot_diagnostics(args...; opts...)
        error(isdefined(Main, :CairoMakie) ? "Invalid method call" : "Import CairoMakie to enable plotting")
    end
    include("data_loading.jl")
    include("preprocessing.jl")
    include("lineage_processing.jl")
end
