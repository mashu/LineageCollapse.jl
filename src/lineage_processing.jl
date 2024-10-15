# Abstract types for the abstraction
abstract type DistanceMetric end
abstract type NormalizedDistanceMetric end
abstract type ClusteringMethod end

# Concrete types for distance metrics
struct HammingDistance <: DistanceMetric end
struct NormalizedHammingDistance <: NormalizedDistanceMetric end
struct LevenshteinDistance <: DistanceMetric end

# Concrete type for hierarchical clustering
struct HierarchicalClustering <: ClusteringMethod
    cutoff::Float64
end

"""
    compute_distance(metric::DistanceMetric, x::LongDNA{4}, y::LongDNA{4})::Float64

Compute the distance between two `LongDNA{4}` sequences using the specified distance metric.
"""
function compute_distance(::HammingDistance, x::LongDNA{4}, y::LongDNA{4})::Float64
    @assert length(x) == length(y)
#    return sum(count_ones.(x.data .⊻ y.data)) / 2
    return evaluate(Hamming(), String(x), String(y))
end

"""
    compute_distance(metric::DistanceMetric, x::LongDNA{4}, y::LongDNA{4})::Float64

Compute the distance between two `LongDNA{4}` sequences using the specified distance metric.
"""
function compute_distance(::NormalizedHammingDistance, x::LongDNA{4}, y::LongDNA{4})::Float64
    @assert length(x) == length(y)
#    return sum(count_ones.(x.data .⊻ y.data)) / 2 / max(length(x), length(y))
    return evaluate(Hamming(), String(x), String(y)) / max(length(x), length(y))
end

"""
    compute_distance(metric::DistanceMetric, x::LongDNA{4}, y::LongDNA{4})::Float64

Compute the distance between two `LongDNA{4}` sequences using the specified distance metric.
"""
function compute_distance(::LevenshteinDistance, x::LongDNA{4}, y::LongDNA{4})::Float64
    return evaluate(Levenshtein(), String(x), String(y))
end

"""
    compute_pairwise_distance(metric::Union{DistanceMetric, NormalizedDistanceMetric}, sequences::Vector{LongDNA{4}})::Matrix{Float64}

Compute pairwise distances between sequences using the specified distance metric.
"""
function compute_pairwise_distance(metric::Union{DistanceMetric, NormalizedDistanceMetric}, sequences::Vector{LongDNA{4}})::Matrix{Float64}
    n = length(sequences)
    dist_matrix = zeros(Float64, n, n)

    Threads.@threads for i in 1:n
        for j in i+1:n
            dist = compute_distance(metric, sequences[i], sequences[j])
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
        end
    end
    return dist_matrix
end

"""
    perform_clustering(method::HierarchicalClustering, dist_matrix::Matrix{Float64})::Vector{Int}

Perform hierarchical clustering on the input distance matrix using the specified method.
"""
function perform_clustering(method::HierarchicalClustering, linkage::Symbol, dist_matrix::Matrix{Float64})::Vector{Int}
    hclusters = hclust(dist_matrix, linkage=linkage)
    return cutree(hclusters, h=method.cutoff)
end

"""
    process_lineages(df::DataFrame; 
                     distance_metric::DistanceMetric = HammingDistance(),
                     clustering_method::ClusteringMethod = HierarchicalClustering(0.1),
                     cdr3_ratio::Float64 = 0.5)::DataFrame

Process lineages in the input DataFrame using specified distance metric and clustering method.

# Arguments
- `df::DataFrame`: Input DataFrame.
- `distance_metric::DistanceMetric`: Distance metric to use for sequence comparison.
- `clustering_method::ClusteringMethod`: Clustering method to use.
- `cdr3_ratio::Float64`: Minimum ratio of CDR3 frequency to consider.

# Returns
- `DataFrame`: Processed DataFrame with lineage information.
"""
function process_lineages(df::DataFrame; 
                          distance_metric::Union{DistanceMetric, NormalizedDistanceMetric} = NormalizedHammingDistance(),
                          clustering_method::ClusteringMethod = HierarchicalClustering(0.1),
                          cdr3_ratio::Float64 = 0.0,
                          linkage::Symbol = :single)::DataFrame
    grouped = groupby(df, [:v_call_first, :j_call_first, :cdr3_length])
    processed_groups = Vector{DataFrame}()

    prog = Progress(length(grouped), desc="Processing lineages")
    for (group_id, group) in enumerate(grouped)
        next!(prog)
        if nrow(group) > 1
            sequences = LongDNA{4}.(group.cdr3)
            dist_matrix = compute_pairwise_distance(distance_metric, sequences)
            group[!, :cluster] = perform_clustering(clustering_method, linkage, dist_matrix)
        else
            group[!, :cluster] .= 1
        end
        
        group[!, :group_id] .= group_id

        cluster_grouped = groupby(group, :cluster)
        for cgroup in cluster_grouped
            cgroup[!, :cluster_size] .= nrow(cgroup)
            cgroup = transform(groupby(cgroup, [:v_call_first, :j_call_first, :cluster, :cdr3_length, :cdr3, :d_region, :cluster_size, :group_id]), nrow => :cdr3_count)
            transform!(groupby(cgroup, :cluster), :cdr3_count => maximum => :max_cdr3_count)
            transform!(groupby(cgroup, :cluster), [:cdr3_count, :max_cdr3_count] => ((count, max_count) -> count ./ max_count) => :cdr3_frequency)
            filter!(row -> row.cdr3_frequency >= cdr3_ratio, cgroup)
            push!(processed_groups, cgroup)
        end
    end
    finish!(prog)

    result = vcat(processed_groups...)
    result[!, :lineage_id] = groupindices(groupby(result, [:group_id, :cluster]))

    return result
end