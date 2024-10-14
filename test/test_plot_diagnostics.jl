using Test
using LineageCollapse

@testset "Plot Diagnostics" begin
    @testset "CairoMakie not imported" begin
        # Ensure CairoMakie is not imported
        if isdefined(Main, :CairoMakie)
            @warn "CairoMakie is already imported, which may affect the test results"
        end

        # Test that the function throws an error with the correct message
        @test_throws ErrorException("Import CairoMakie to enable plotting") plot_diagnostics()
        
        # Test with arguments
        @test_throws ErrorException("Import CairoMakie to enable plotting") plot_diagnostics(rand(10))
        
        # Test with keyword arguments
        @test_throws ErrorException("Import CairoMakie to enable plotting") plot_diagnostics(color="red")
    end

    @testset "CairoMakie imported" begin
        # Create a temporary module to simulate CairoMakie being imported
        @eval Main module TestPlotDiagnostics
            using LineageCollapse
            
            # Simulate CairoMakie being imported
            LineageCollapse.eval(:(CairoMakie = 1))

            using Test
            @test_throws ErrorException("Invalid method call") plot_diagnostics()
            @test_throws ErrorException("Invalid method call") plot_diagnostics(rand(10))
            @test_throws ErrorException("Invalid method call") plot_diagnostics(color="red")
        end
    end
end