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

"""
    process_lineages(df::DataFrame;
                    distance_metric::Union{DistanceMetric, NormalizedDistanceMetric} = NormalizedHammingDistance(),
                    clustering_method::ClusteringMethod = HierarchicalClustering(0.1),
                    linkage::Symbol = :single)::DataFrame

Process lineages from a DataFrame of CDR3 sequences.
"""
function process_lineages(df::DataFrame;
                          distance_metric::Union{DistanceMetric, NormalizedDistanceMetric} = NormalizedHammingDistance(),
                          clustering_method::ClusteringMethod = HierarchicalClustering(0.2),
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
    collapse_lineages(df::DataFrame, clone_frequency_threshold::Float64, collapse_strategy::Symbol=:hardest)

Collapse lineages in a DataFrame based on clone frequency and a specified collapse strategy.

# Arguments
- `df::DataFrame`: Input DataFrame containing lineage data. Must have columns [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3].
- `clone_frequency_threshold::Float64`: Minimum frequency threshold for clones (0.0 to 1.0).
- `collapse_strategy::Symbol=:hardest`: Strategy for collapsing lineages. Options are:
  - `:hardest`: Select only the most frequent clone for each lineage.
  - `:soft`: Select all clones that meet or exceed the `clone_frequency_threshold`.

# Returns
- `DataFrame`: Collapsed lineage data containing a new column:
  - `clone_frequency`: Represents the relative frequency of each clone within its lineage,
    calculated as (count of specific clone) / (total sequences in lineage).
    A clone is defined by unique combination of D region, J call, V call, and CDR3 sequence.
    Values range from 0.0 to 1.0, with higher values indicating more abundant clones in the lineage.

# Example
```julia
lineages = DataFrame(...)  # Your input data
collapsed = collapse_lineages(lineages, 0.1, :soft)
```
"""
function collapse_lineages(df::DataFrame, clone_frequency_threshold::Float64, collapse_strategy::Symbol=:hardest)
    if !(0.0 <= clone_frequency_threshold <= 1.0)
        throw(ArgumentError("clone_frequency_threshold must be between 0.0 and 1.0"))
    end
    if !(collapse_strategy in [:hardest, :soft])
        throw(ArgumentError("Invalid collapse strategy. Use :hardest or :soft."))
    end

    grouped = groupby(df, [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3])
    counted = combine(grouped, nrow => :sequence_count)

    lineage_grouped = groupby(counted, :lineage_id)
    with_frequency = transform(lineage_grouped, :sequence_count => (x -> x ./ sum(x)) => :clone_frequency)

    df_with_freq = leftjoin(df, with_frequency,
        on=[:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3],
        makeunique=true)

    collapsed = if collapse_strategy == :hardest
        combine(groupby(df_with_freq, :lineage_id)) do group
            row_idx = argmax(group.clone_frequency)  # Single highest frequency
            group[row_idx:row_idx, :]
        end
    else
        combine(groupby(df_with_freq, :lineage_id)) do group
            group[group.clone_frequency .>= clone_frequency_threshold, :]
        end
    end

    return collapsed
end