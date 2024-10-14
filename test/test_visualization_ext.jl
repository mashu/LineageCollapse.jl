using Test
using DataFrames
using LineageCollapse
using CairoMakie

@testset "VisualizationExt" begin
    @testset "plot_diagnostics function" begin
        # Helper function to check figure content
        function check_figure_content(fig)
            @test fig isa CairoMakie.Figure
            @test length(fig.content) == 9  # 4 subplots + 1 title + 4 subplot labels
            @test any(x -> x isa Axis, fig.content)
            @test any(x -> x isa Label, fig.content)
        end

        # Test with all required columns
        df_complete = DataFrame(
            cluster_size = rand(1:100, 1000),
            cdr3_length = rand(5:30, 1000),
            cdr3_frequency = rand(1000),
            v_call_first = rand(["IGHV1-1", "IGHV1-2", "IGHV1-3", "IGHV2-1", "IGHV3-1"], 1000)
        )
        
        fig_complete = plot_diagnostics(df_complete)
        check_figure_content(fig_complete)

        # Test with missing cluster_size
        df_no_cluster = select(df_complete, Not(:cluster_size))
        fig_no_cluster = plot_diagnostics(df_no_cluster)
        check_figure_content(fig_no_cluster)

        # Test with missing cdr3_length
        df_no_cdr3_length = select(df_complete, Not(:cdr3_length))
        fig_no_cdr3_length = plot_diagnostics(df_no_cdr3_length)
        check_figure_content(fig_no_cdr3_length)

        # Test with missing cdr3_frequency
        df_no_cdr3_freq = select(df_complete, Not(:cdr3_frequency))
        fig_no_cdr3_freq = plot_diagnostics(df_no_cdr3_freq)
        check_figure_content(fig_no_cdr3_freq)

        # Test with missing v_call_first
        df_no_v_call = select(df_complete, Not(:v_call_first))
        fig_no_v_call = plot_diagnostics(df_no_v_call)
        check_figure_content(fig_no_v_call)

        # Test with empty DataFrame
        df_empty = DataFrame()
        fig_empty = plot_diagnostics(df_empty)
        check_figure_content(fig_empty)

        # Test with DataFrame containing only some columns
        df_partial = DataFrame(
            cluster_size = rand(1:100, 100),
            v_call_first = rand(["IGHV1-1", "IGHV1-2", "IGHV1-3"], 100)
        )
        fig_partial = plot_diagnostics(df_partial)
        check_figure_content(fig_partial)
    end
end