# Abstract types for the abstraction
abstract type DistanceMetric end
abstract type NormalizedDistanceMetric end
abstract type ClusteringMethod end

# Concrete types for distance metrics
struct HammingDistance <: DistanceMetric end
struct NormalizedHammingDistance <: NormalizedDistanceMetric end
struct LevenshteinDistance <: DistanceMetric end

"""
    HierarchicalClustering(cutoff::Float32)

A type representing hierarchical clustering with a cutoff.

# Arguments
- `cutoff::Float32`: The cutoff value for the clustering, below which clusters are merged. Higher values result in fewer clusters.
"""
struct HierarchicalClustering <: ClusteringMethod
    cutoff::Float32
end

"""
    compute_distance(metric::DistanceMetric, x::LongDNA{4}, y::LongDNA{4})::Float32

Compute the distance between two `LongDNA{4}` sequences using the specified distance metric.
"""
function compute_distance(::HammingDistance, x::LongDNA{4}, y::LongDNA{4})::Float32
    @assert length(x) == length(y)
    return mismatches(x, y)
end

"""
    compute_distance(metric::DistanceMetric, x::LongDNA{4}, y::LongDNA{4})::Float32

Compute the distance between two `LongDNA{4}` sequences using the specified distance metric.
"""
function compute_distance(::NormalizedHammingDistance, x::LongDNA{4}, y::LongDNA{4})::Float32
    @assert length(x) == length(y)
    return mismatches(x, y) / length(x)
end

"""
    compute_distance(metric::DistanceMetric, x::LongDNA{4}, y::LongDNA{4})::Float32

Compute the distance between two `LongDNA{4}` sequences using the specified distance metric.
"""
function compute_distance(::LevenshteinDistance, x::LongDNA{4}, y::LongDNA{4})::Float32
    return evaluate(Levenshtein(), String(x), String(y))
end

"""
    compute_pairwise_distance(metric::Union{DistanceMetric, NormalizedDistanceMetric}, sequences::Vector{LongDNA{4}})::Matrix{Float32}

Compute pairwise distances between sequences using the specified distance metric.
"""
function compute_pairwise_distance(
    metric::M, 
    sequences::AbstractVector{S}
)::Matrix{Float32} where {M <: Union{DistanceMetric, NormalizedDistanceMetric}, S <: LongSequence{DNAAlphabet{4}}}
    n = length(sequences)
    @assert n > 0 "No sequences provided for distance calculation"
    dist_matrix = spzeros(Float32, n, n)

    Threads.@threads for i in 1:n
        for j in i+1:n
            @inbounds dist = compute_distance(metric, sequences[i], sequences[j])
            @inbounds dist_matrix[i, j] = dist
        end
    end
    return Symmetric(dist_matrix)
end

"""
    perform_clustering(method::HierarchicalClustering, linkage::Symbol, dist_matrix::T)::Vector{Int} where T <: AbstractMatrix

Perform hierarchical clustering on the distance matrix using the specified method and linkage.
"""
function perform_clustering(method::HierarchicalClustering, linkage::Symbol, dist_matrix::T)::Vector{Int} where T <: AbstractMatrix
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
    # Convert upfront to LongDNA{4} for performance
    df.cdr3 = LongDNA{4}.(df.cdr3)

    grouped = groupby(df, [:v_call_first, :j_call_first, :cdr3_length])
    processed_groups = Vector{DataFrame}()

    prog = Progress(length(grouped), desc="Processing lineages")
    @inbounds for (group_id, group) in enumerate(grouped)
        next!(prog)
        if nrow(group) > 1
            dist_matrix = compute_pairwise_distance(distance_metric, group.cdr3)
            group[!, :cluster] = perform_clustering(clustering_method, linkage, dist_matrix)
        else
            group[!, :cluster] .= 1
        end

        group[!, :group_id] .= group_id

        cluster_grouped = groupby(group, :cluster)
        @inbounds for cgroup in cluster_grouped
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

"""
    collapse_lineages(df::DataFrame, cdr3_frequency_threshold::Float64, collapse_strategy::Symbol=:hardest)

Collapse lineages in a DataFrame based on CDR3 sequence frequency and a specified collapse strategy.

# Arguments
- `df::DataFrame`: Input DataFrame containing lineage data. Must have columns [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3].
- `cdr3_frequency_threshold::Float64`: Minimum frequency threshold for CDR3 sequences (0.0 to 1.0).
- `collapse_strategy::Symbol=:hardest`: Strategy for collapsing lineages. Options are:
  - `:hardest`: Select only the most frequent sequence for each lineage.
  - `:soft`: Select all sequences that meet or exceed the `cdr3_frequency_threshold`.

# Returns
- `DataFrame`: Collapsed lineage data.

# Example
```julia
lineages = DataFrame(...)  # Your input data
collapsed = collapse_lineages(lineages, 0.1, :soft)
```
"""
function collapse_lineages(df::DataFrame, cdr3_frequency_threshold::Float64, collapse_strategy::Symbol=:hardest)
    if !(0.0 <= cdr3_frequency_threshold <= 1.0)
        throw(ArgumentError("cdr3_frequency_threshold must be between 0.0 and 1.0"))
    end
    if !(collapse_strategy in [:hardest, :soft])
        throw(ArgumentError("Invalid collapse strategy. Use :hardest or :soft."))
    end

    # Group by the specified columns
    grouped = groupby(df, [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3])

    # Count occurrences of each unique combination
    counted = combine(grouped, nrow => :count)

    # Calculate frequency within each lineage
    lineage_grouped = groupby(counted, :lineage_id)
    with_frequency = transform(lineage_grouped, :count => (x -> x ./ sum(x)) => :frequency)

    # Filter based on the threshold and collapse strategy
    function filter_lineage(group)
        if collapse_strategy == :hardest
            # Pick the single most frequent sequence
            return group[argmax(group.frequency), :]
        else
            # Pick all sequences above the threshold
            return group[group.frequency .>= cdr3_frequency_threshold, :]
        end
    end

    collapsed = combine(groupby(with_frequency, :lineage_id), filter_lineage)

    return collapsed
end