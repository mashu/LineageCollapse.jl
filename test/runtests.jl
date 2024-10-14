using Test
using LineageCollapse

@testset "LineageCollapse.jl" begin
    include("test_data_loading.jl")
    include("test_preprocessing.jl")
    include("test_lineage_processing.jl")
    include("test_plot_diagnostics.jl")
end