# Abstract types for the abstraction
abstract type DistanceMetric end
abstract type NormalizedDistanceMetric end
abstract type ClusteringMethod end

# Concrete types for distance metrics
struct HammingDistance <: DistanceMetric end
struct NormalizedHammingDistance <: NormalizedDistanceMetric end
struct LevenshteinDistance <: DistanceMetric end

"""
    HierarchicalClustering(cutoff::Float64)

A type representing hierarchical clustering with a cutoff.

# Arguments
- `cutoff::Float64`: The cutoff value for the clustering, below which clusters are merged. Higher values result in fewer clusters.
"""
struct HierarchicalClustering <: ClusteringMethod
    cutoff::Float64
end

"""
    compute_distance(metric::DistanceMetric, x::LongDNA{4}, y::LongDNA{4})::Float64

Compute the distance between two `LongDNA{4}` sequences using the specified distance metric.
"""
function compute_distance(::HammingDistance, x::LongDNA{4}, y::LongDNA{4})::Float64
    @assert length(x) == length(y)
    return evaluate(Hamming(), String(x), String(y))
end

"""
    compute_distance(metric::DistanceMetric, x::LongDNA{4}, y::LongDNA{4})::Float64

Compute the distance between two `LongDNA{4}` sequences using the specified distance metric.
"""
function compute_distance(::NormalizedHammingDistance, x::LongDNA{4}, y::LongDNA{4})::Float64
    @assert length(x) == length(y)
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
    perform_clustering(method::HierarchicalClustering, linkage::Symbol, dist_matrix::Matrix{Float64})::Vector{Int}

Perform hierarchical clustering on the distance matrix using the specified method and linkage.
"""
function perform_clustering(method::HierarchicalClustering, linkage::Symbol, dist_matrix::Matrix{Float64})::Vector{Int}
    hclusters = hclust(dist_matrix, linkage=linkage)
    return cutree(hclusters, h=method.cutoff)
end

"""
    process_lineages(df::DataFrame; 
                    distance_metric::Union{DistanceMetric, NormalizedDistanceMetric} = NormalizedHammingDistance(),
                    clustering_method::ClusteringMethod = HierarchicalClustering(0.1),
                    linkage::Symbol = :single)::DataFrame

Process lineages from a DataFrame of CDR3 sequences.
"""
function process_lineages(df::DataFrame; 
                          distance_metric::Union{DistanceMetric, NormalizedDistanceMetric} = NormalizedHammingDistance(),
                          clustering_method::ClusteringMethod = HierarchicalClustering(0.1),
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
            push!(processed_groups, cgroup)
        end
    end
    finish!(prog)

    result = vcat(processed_groups...)
    result[!, :lineage_id] = groupindices(groupby(result, [:group_id, :cluster]))

    return result
end