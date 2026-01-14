# Abstract types for the abstraction
abstract type DistanceMetric end
abstract type NormalizedDistanceMetric end
abstract type ClusteringMethod end
abstract type CollapseStrategy end

# Concrete types for distance metrics
struct HammingDistance <: DistanceMetric end
struct NormalizedHammingDistance <: NormalizedDistanceMetric end
struct LevenshteinDistance <: DistanceMetric end

# Collapse strategies
struct Hardest <: CollapseStrategy end
struct Soft{T <: AbstractFloat} <: CollapseStrategy
    cutoff::T
    function Soft{T}(cutoff::T) where {T <: AbstractFloat}
        if !(0.0 <= cutoff <= 1.0)
            throw(ArgumentError("clone_frequency_threshold must be between 0.0 and 1.0"))
        end
        return new{T}(cutoff)
    end
end

Soft(cutoff::AbstractFloat) = Soft{typeof(cutoff)}(cutoff)
Soft(cutoff::Integer) = Soft(float(cutoff))

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
    dist_matrix = zeros(Float32, n, n)

    Threads.@threads for i in 1:n
        for j in i+1:n
            @inbounds dist = compute_distance(metric, sequences[i], sequences[j])
            @inbounds dist_matrix[i, j] = dist
        end
    end
    return dist_matrix
end

"""
    perform_clustering(method::HierarchicalClustering, linkage::Symbol, dist_matrix::T)::Vector{Int} where T <: AbstractMatrix

Perform hierarchical clustering on the distance matrix using the specified method and linkage.
"""
function perform_clustering(method::HierarchicalClustering, linkage::Symbol, dist_matrix::T)::Vector{Int} where T <: AbstractMatrix
    hclusters = hclust(dist_matrix, linkage=linkage, uplo=:U)
    return cutree(hclusters, h=method.cutoff)
end

_cdr3_threshold_params(threshold::Integer) = begin
    if threshold < 0
        throw(ArgumentError("cdr3_mismatch_threshold must be non-negative"))
    end
    return HammingDistance(), HierarchicalClustering(Float32(threshold))
end

_cdr3_threshold_params(threshold::AbstractFloat) = begin
    if !(0.0 <= threshold <= 1.0)
        throw(ArgumentError("cdr3_mismatch_threshold as a fraction must be between 0.0 and 1.0"))
    end
    return NormalizedHammingDistance(), HierarchicalClustering(Float32(threshold))
end

"""
    process_lineages(df::DataFrame, cdr3_mismatch_threshold::Union{Integer,AbstractFloat};
                     linkage::Symbol = :single)::DataFrame

Process lineages from a DataFrame of CDR3 sequences using a mismatch threshold.

If `cdr3_mismatch_threshold` is an `Integer`, it is interpreted as an absolute
number of allowed mismatches. If it is a floating-point value, it is interpreted
as a fraction of the CDR3 length (0.0 to 1.0). Use an integer literal (e.g. `1`)
to request absolute mismatches when `1.0` would be ambiguous.
"""
function process_lineages(df::DataFrame, cdr3_mismatch_threshold::Union{Integer,AbstractFloat};
                          linkage::Symbol = :single)::DataFrame
    distance_metric, clustering_method = _cdr3_threshold_params(cdr3_mismatch_threshold)
    return process_lineages(df;
        distance_metric=distance_metric,
        clustering_method=clustering_method,
        linkage=linkage,
    )
end

"""
    process_lineages(df::DataFrame;
                    distance_metric::Union{DistanceMetric, NormalizedDistanceMetric} = HammingDistance(),
                    clustering_method::ClusteringMethod = HierarchicalClustering(1.0),
                    linkage::Symbol = :single)::DataFrame

Process lineages from a DataFrame of CDR3 sequences.
"""
function process_lineages(df::DataFrame;
                          distance_metric::Union{DistanceMetric, NormalizedDistanceMetric} = HammingDistance(),
                          clustering_method::ClusteringMethod = HierarchicalClustering(1.0),
                          linkage::Symbol = :single)::DataFrame
    # Group by VJ combination and CDR3 length first
    groups = groupby(df, [:v_call_first, :j_call_first, :cdr3_length])
    processed_groups = Vector{DataFrame}()

    for group in groups
        # Within each VJ+length group, get unique CDR3s
        unique_dna_data = unique(group, :cdr3)

        if nrow(unique_dna_data) > 1
            dna_seqs = LongDNA{4}.(unique_dna_data.cdr3)
            dist_matrix = compute_pairwise_distance(distance_metric, dna_seqs)

            # Calculate minimum distances for each sequence, excluding zeros on diagonal
            min_distances = zeros(Float32, size(dist_matrix, 1))
            for i in 1:size(dist_matrix, 1)
                non_zero_distances = filter(x -> x > 0, dist_matrix[i,:])
                min_distances[i] = isempty(non_zero_distances) ? 0.0f0 : minimum(non_zero_distances)
            end
            unique_dna_data[!, :min_distance] = min_distances

            unique_dna_data[!, :cluster] = perform_clustering(clustering_method, linkage, dist_matrix)
        else
            unique_dna_data[!, :cluster] .= 1
            unique_dna_data[!, :min_distance] .= 0.0f0  # Single sequence has no distances to others
        end

        # Map clusters and min distances back to all sequences
        group_result = leftjoin(group,
                              select(unique_dna_data, :cdr3, :cluster, :min_distance),
                              on = :cdr3)

        # Process statistics
        transform!(group_result, :cluster => (x -> length.(x)) => :cluster_size)
        group_result = transform(groupby(group_result, [:cluster, :cdr3]), nrow => :cdr3_count)
        transform!(groupby(group_result, :cluster), :cdr3_count => maximum => :max_cdr3_count)
        transform!(groupby(group_result, :cluster),
                  [:cdr3_count, :max_cdr3_count] =>
                  ((count, max_count) -> count ./ max_count) => :cdr3_frequency)

        push!(processed_groups, group_result)
    end

    result = vcat(processed_groups...)
    result[!, :lineage_id] = groupindices(groupby(result, [:v_call_first, :j_call_first, :cdr3_length, :cluster]))

    return result
end

"""
    collapse_lineages(df::DataFrame, strategy::CollapseStrategy=Hardest())

Collapse lineages in a DataFrame based on clone frequency and a specified collapse strategy.

# Arguments
- `df::DataFrame`: Input DataFrame containing lineage data. Must have columns [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3].
- `strategy::CollapseStrategy=Hardest()`: Strategy for collapsing lineages.
  - `Hardest()` selects only the most frequent clone per lineage.
  - `Soft(cutoff)` keeps clones whose frequency within a lineage is at or above `cutoff`.

# Returns
- `DataFrame`: Collapsed lineage data containing a new column:
  - `clone_frequency`: Represents the relative frequency of each clone within its lineage,
    calculated as (count of specific clone) / (total sequences in lineage).
    A clone is defined by unique combination of D region, J call, V call, and CDR3 sequence.
    Values range from 0.0 to 1.0, with higher values indicating more abundant clones in the lineage.

# Example
```julia
lineages = DataFrame(...)  # Your input data
collapsed = collapse_lineages(lineages, Soft(0.1))
```
"""
function collapse_lineages(df::DataFrame, strategy::CollapseStrategy=Hardest())
    grouped = groupby(df, [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3])
    counted = combine(grouped, nrow => :sequence_count)

    lineage_grouped = groupby(counted, :lineage_id)
    with_frequency = transform(lineage_grouped, :sequence_count => (x -> x ./ sum(x)) => :clone_frequency)

    df_with_freq = leftjoin(df, with_frequency,
        on=[:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3],
        makeunique=true)

    collapsed = _collapse_by_strategy(df_with_freq, strategy)

    return collapsed
end

function _collapse_by_strategy(df_with_freq::DataFrame, ::Hardest)
    return combine(groupby(df_with_freq, :lineage_id)) do group
        row_idx = argmax(group.clone_frequency)  # Single highest frequency
        group[row_idx:row_idx, :]
    end
end

function _collapse_by_strategy(df_with_freq::DataFrame, strategy::Soft)
    return combine(groupby(df_with_freq, :lineage_id)) do group
        group[group.clone_frequency .>= strategy.cutoff, :]
    end
end

# Backward-compatible API