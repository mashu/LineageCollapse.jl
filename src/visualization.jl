"""
    plot_diagnostics(df::DataFrame)

Generate diagnostic plots for the processed data.

# Arguments
- `df::DataFrame`: Processed DataFrame with lineage information.

# Returns
- `Plots.Plot`: A composite plot with multiple diagnostic visualizations.
"""
function plot_diagnostics(df::DataFrame)
    p1 = histogram(df.cluster_size, title="Cluster Size Distribution", xlabel="Cluster Size", ylabel="Frequency")
    
    p2 = scatter(df.cdr3_length, df.cluster_size, title="CDR3 Length vs Cluster Size", xlabel="CDR3 Length", ylabel="Cluster Size")
    
    p3 = histogram(df.cdr3_frequency, title="CDR3 Frequency Distribution", xlabel="CDR3 Frequency", ylabel="Count")
    
    v_gene_counts = sort(combine(groupby(df, :v_call_first), nrow => :count), :count, rev=true)[1:10, :]
    p4 = bar(v_gene_counts.v_call_first, v_gene_counts.count, title="Top 10 V Genes", xrotation=90, xlabel="V Gene", ylabel="Count")
    
    return plot(p1, p2, p3, p4, layout=(2,2), size=(1200,800))
end