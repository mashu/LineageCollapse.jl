"""
    hamming(x::LongDNA{4}, y::LongDNA{4})::Float64

Calculate the Hamming distance between two `LongDNA{4}` sequences.

# Returns
- `Float64`: Hamming distance.
"""
function hamming(x::LongDNA{4}, y::LongDNA{4})::Float64
    xor_result = x.data .⊻ y.data
    return sum(count_ones.(xor_result)) / 2
end

"""
    pairwise_hamming(sequences::Vector{LongDNA{4}})::Array{Float64, 2}

Compute the pairwise Hamming distance matrix for a vector of `LongDNA{4}` sequences.

# Returns
- `Array{Float64, 2}`: Symmetric distance matrix.
"""
function pairwise_hamming(sequences::Vector{LongDNA{4}})
    n = length(sequences)
    dist_matrix = zeros(Float64, n, n)

    Threads.@threads for i in 1:n
        for j in i+1:n
            dist = hamming(sequences[i], sequences[j])
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
        end
    end

    return dist_matrix
end

"""
    process_lineages(df::DataFrame; 
                     cutoff_ratio::Float64=0.1, 
                     allele_ratio::Float64=0.5)::DataFrame

Process lineages in the input DataFrame.

# Arguments
- `df::DataFrame`: Input DataFrame.
- `cutoff_ratio::Float64=0.1`: Ratio to determine the cutoff height for hierarchical clustering.
- `allele_ratio::Float64=0.5`: Minimum ratio of CDR3 frequency to consider.

# Returns
- `DataFrame`: Processed DataFrame with lineage information.
"""
function process_lineages(df::DataFrame; 
                          cutoff_ratio::Float64=0.1, 
                          allele_ratio::Float64=0.5)::DataFrame
    grouped = groupby(df, [:v_call_first, :j_call_first, :cdr3_length])
    processed_groups = Vector{DataFrame}()

    @showprogress "Processing lineages" for group in grouped
        if nrow(group) > 1
            dist = pairwise_hamming(LongDNA{4}.(group.cdr3))
            hclusters = hclust(dist, linkage=:average)
            maximum_height = maximum(hclusters.heights)
            cutoff = maximum_height * cutoff_ratio
            group[!, :cluster] = cutree(hclusters, h=cutoff)
        else
            group[!, :cluster] .= 1
        end

        cluster_grouped = groupby(group, :cluster)
        for cgroup in cluster_grouped
            cgroup[!, :cluster_size] .= nrow(cgroup)
            cgroup = combine(groupby(cgroup, [:v_call_first, :j_call_first, :cluster, :cdr3_length, :cdr3, :d_region, :cluster_size]), nrow => :cdr3_count)
            transform!(groupby(cgroup, :cluster), :cdr3_count => maximum => :max_cdr3_count)
            transform!(groupby(cgroup, :cluster), [:cdr3_count, :max_cdr3_count] => ((count, max_count) -> count ./ max_count) => :cdr3_frequency)
            filter!(row -> row.cdr3_frequency > allele_ratio, cgroup)
            push!(processed_groups, cgroup)
        end
    end

    return vcat(processed_groups...)
end
