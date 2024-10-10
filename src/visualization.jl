"""
    plot_diagnostics(df::DataFrame)

Generate professional diagnostic plots for the processed data using CairoMakie.

# Arguments
- `df::DataFrame`: Processed DataFrame with lineage information.

# Returns
- `CairoMakie.Figure`: A composite figure with multiple diagnostic visualizations.
"""
function plot_diagnostics(df::DataFrame)
    fig = Figure(size=(1200, 1000), fontsize=12)

    # Cluster Size Distribution
    if hasproperty(df, :cluster_size)
        ax1 = Axis(fig[1, 1], title="Cluster Size Distribution", xlabel="Cluster Size", ylabel="Frequency")
        hist!(ax1, df.cluster_size, bins=50, color=:skyblue, strokecolor=:white, strokewidth=1)
    else
        ax1 = Axis(fig[1, 1], title="Cluster Size Distribution Not Available")
    end

    # CDR3 Length vs Cluster Size
    if hasproperty(df, :cdr3_length) && hasproperty(df, :cluster_size)
        ax2 = Axis(fig[1, 2], title="CDR3 Length vs Cluster Size", xlabel="CDR3 Length", ylabel="Cluster Size")
        scatter!(ax2, df.cdr3_length, df.cluster_size, color=:darkblue, markersize=4, alpha=0.5)
    else
        ax2 = Axis(fig[1, 2], title="CDR3 Length vs Cluster Size Not Available")
    end

    # CDR3 Frequency Distribution
    if hasproperty(df, :cdr3_frequency)
        ax3 = Axis(fig[2, 1], title="CDR3 Frequency Distribution", xlabel="CDR3 Frequency", ylabel="Count")
        hist!(ax3, df.cdr3_frequency, bins=50, color=:lightgreen, strokecolor=:white, strokewidth=1)
    else
        ax3 = Axis(fig[2, 1], title="CDR3 Frequency Distribution Not Available")
    end

    # Top 10 V Genes
    ax4 = Axis(fig[2, 2], title="Top 10 V Genes", xlabel="V Gene", ylabel="Count")
    if hasproperty(df, :v_call_first)
        try
            v_gene_counts = sort(combine(groupby(df, :v_call_first), nrow => :count), :count, rev=true)
            if nrow(v_gene_counts) > 10
                v_gene_counts = v_gene_counts[1:10, :]
            end
            barplot!(ax4, v_gene_counts.count, color=:orange)
            ax4.xticks = (1:nrow(v_gene_counts), v_gene_counts.v_call_first)
            ax4.xticklabelrotation = π/3
        catch e
            @warn "Error processing V gene counts: $e"
            ax4.title = "Top 10 V Genes (Error in Processing)"
        end
    else
        ax4.title = "Top 10 V Genes Not Available"
    end

    # Rotate x-axis labels for all plots
    for ax in [ax1, ax2, ax3, ax4]
        ax.xticklabelrotation = π/4
    end

    # Add a title to the entire figure
    Label(fig[0, :], "Diagnostic Plots for Lineage Collapse", fontsize=20)

    # Adjust layout
    for (label, layout) in zip(["A", "B", "C", "D"], [fig[1,1], fig[1,2], fig[2,1], fig[2,2]])
        Label(layout[1, 1, TopLeft()], label,
              fontsize = 26,
              padding = (6, 6, 6, 6),
              halign = :right)
    end

    return fig
end