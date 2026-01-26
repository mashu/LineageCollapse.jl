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

#=
    Lineage Aggregation

    Multiple dispatch handles whether to aggregate based on strategy type.
    Hardest() aggregates count and computes nVDJ_nt.
    Soft() does not aggregate.
=#

struct LineageAggregates
    data::DataFrame
end

struct NoAggregates end

"""
    compute_aggregates(::Hardest, df::AbstractDataFrame) -> LineageAggregates

Compute aggregated statistics for each lineage when using Hardest strategy.
"""
function compute_aggregates(::Hardest, df::AbstractDataFrame)
    has_count = :count ∈ propertynames(df)
    has_vdj_nt = :vdj_nt ∈ propertynames(df)
    
    agg = combine(groupby(df, :lineage_id)) do group
        total_count = has_count ? sum(skipmissing(group.count)) : nrow(group)
        n_vdj_nt = has_vdj_nt ? length(unique(skipmissing(group.vdj_nt))) : missing
        (lineage_count_sum = total_count, nVDJ_nt = n_vdj_nt)
    end
    return LineageAggregates(agg)
end

"""
    compute_aggregates(::Soft, ::AbstractDataFrame) -> NoAggregates

Soft strategy does not aggregate - returns sentinel type.
"""
compute_aggregates(::Soft, ::AbstractDataFrame) = NoAggregates()

"""
    apply_aggregates(collapsed::DataFrame, agg::LineageAggregates) -> DataFrame

Apply precomputed aggregates to the collapsed DataFrame.
"""
function apply_aggregates(collapsed::DataFrame, agg::LineageAggregates)
    result = leftjoin(collapsed, agg.data, on=:lineage_id)
    
    if :count ∈ propertynames(result) && :lineage_count_sum ∈ propertynames(result)
        result.count = result.lineage_count_sum
        select!(result, Not(:lineage_count_sum))
    elseif :lineage_count_sum ∈ propertynames(result)
        rename!(result, :lineage_count_sum => :count)
    end
    
    return result
end

"""
    apply_aggregates(collapsed::DataFrame, ::NoAggregates) -> DataFrame

No-op when there are no aggregates to apply.
"""
apply_aggregates(collapsed::DataFrame, ::NoAggregates) = collapsed

#=
    Clone Frequency Computation
=#

"""
    clone_frequency_table(df::DataFrame) -> DataFrame

Compute clone frequency for each unique clone within each lineage.
"""
function clone_frequency_table(df::DataFrame)
    grouped = groupby(df, [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3])
    counted = combine(grouped, nrow => :sequence_count)
    lineage_grouped = groupby(counted, :lineage_id)
    return transform(lineage_grouped, :sequence_count => (x -> x ./ sum(x)) => :clone_frequency)
end

#=
    VDJ Count Computation
=#

"""
    add_vdj_count(df::DataFrame) -> DataFrame

Add vdj_count column showing the count of each unique VDJ sequence within each clone.
"""
function add_vdj_count(df::DataFrame)
    weight_col = :count ∈ propertynames(df) ? :count : nothing
    
    grouped = if weight_col === nothing
        combine(
            groupby(df, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3, :vdj_nt]),
            nrow => :vdj_count,
        )
    else
        combine(
            groupby(df, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3, :vdj_nt]),
            :count => sum => :vdj_count,
        )
    end

    by_clone = combine(
        groupby(grouped, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3]),
        :vdj_count => maximum => :vdj_count,
    )

    return leftjoin(df, by_clone,
        on=[:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3],
        makeunique=true)
end

#=
    Group Collapsing - Multiple Dispatch on Strategy
=#

"""
    collapse_group(::Hardest, group::AbstractDataFrame, tie_breaker::AbstractTieBreaker, tie_atol::Real)

Collapse a lineage group using Hardest strategy - select the most frequent clone.
"""
function collapse_group(::Hardest, group::AbstractDataFrame, tie_breaker::AbstractTieBreaker, tie_atol::Real)
    max_freq = maximum(group.clone_frequency)
    
    candidates = if tie_atol > 0
        group[abs.(group.clone_frequency .- max_freq) .<= tie_atol, :]
    else
        group[group.clone_frequency .== max_freq, :]
    end
    
    return select_representative(tie_breaker, group, candidates)
end

"""
    collapse_group(strategy::Soft, group::AbstractDataFrame, ::AbstractTieBreaker, ::Real)

Collapse a lineage group using Soft strategy - keep clones above frequency cutoff.
"""
function collapse_group(strategy::Soft, group::AbstractDataFrame, ::AbstractTieBreaker, ::Real)
    return group[group.clone_frequency .>= strategy.cutoff, :]
end

#=
    Representative Selection - Multiple Dispatch on TieBreaker
=#

"""
    select_representative(::MostCommonVdjNtTieBreaker, group::AbstractDataFrame, ::AbstractDataFrame)

Select representative using igdiscover's logic: most common VDJ_nt weighted by count.
Uses the entire group (not just frequency-tied candidates) to match igdiscover behavior.
"""
function select_representative(::MostCommonVdjNtTieBreaker, group::AbstractDataFrame, ::AbstractDataFrame)
    if :vdj_nt ∉ propertynames(group)
        throw(ArgumentError("vdj_nt column is required for ByMostCommonVdjNt() tie-breaker."))
    end
    
    n = nrow(group)
    
    # Match igdiscover's behavior: if n <= 2, just take the first row
    # This ensures deterministic behavior when there are ties
    if n <= 2
        return group[1:1, :]
    end
    
    # For n > 2, use most common VDJ_nt weighted by count
    count_col = :count ∈ propertynames(group) ? group.count : ones(Int, nrow(group))
    
    vdj_counts = Dict{String,Int}()
    for (vdj_nt, cnt) in zip(group.vdj_nt, count_col)
        if !ismissing(vdj_nt)
            vdj_counts[vdj_nt] = get(vdj_counts, vdj_nt, 0) + cnt
        end
    end
    
    if isempty(vdj_counts)
        return group[1:1, :]
    end
    
    most_common_vdj_nt = argmax(vdj_counts)
    
    for i in 1:nrow(group)
        if !ismissing(group.vdj_nt[i]) && group.vdj_nt[i] == most_common_vdj_nt
            return group[i:i, :]
        end
    end
    
    return group[1:1, :]
end

"""
    select_representative(tie_breaker::TieBreaker, ::AbstractDataFrame, candidates::AbstractDataFrame)

Select representative from candidates using TieBreaker sorting criteria.
"""
function select_representative(tie_breaker::TieBreaker, ::AbstractDataFrame, candidates::AbstractDataFrame)
    return select_hardest_candidate(candidates, tie_breaker)
end

#=
    VDJ_NT Requirement Check - Multiple Dispatch
=#

"""
    requires_vdj_nt(::MostCommonVdjNtTieBreaker) -> Bool
"""
requires_vdj_nt(::MostCommonVdjNtTieBreaker) = true

"""
    requires_vdj_nt(tie_breaker::TieBreaker) -> Bool
"""
requires_vdj_nt(tie_breaker::TieBreaker) = any(col == :vdj_count for (col, _) in tie_breaker.criteria)

"""
    requires_vdj_count(::MostCommonVdjNtTieBreaker) -> Bool
"""
requires_vdj_count(::MostCommonVdjNtTieBreaker) = false

"""
    requires_vdj_count(tie_breaker::TieBreaker) -> Bool
"""
requires_vdj_count(tie_breaker::TieBreaker) = any(col == :vdj_count for (col, _) in tie_breaker.criteria)

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
    collapse_lineages(df::DataFrame, strategy::CollapseStrategy=Hardest();
                      tie_breaker::AbstractTieBreaker=ByVdjCount(),
                      tie_atol::Real=0.0)

Collapse lineages in a DataFrame based on clone frequency and a specified collapse strategy.

When using `Hardest()` strategy (one representative per lineage), the function automatically:
- Sums the `count` column across all members of each lineage
- Adds `nVDJ_nt` column with the number of unique VDJ_nt sequences per lineage

This matches igdiscover's clonotypes output format.

# Arguments
- `df::DataFrame`: Input DataFrame containing lineage data. Must have columns [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3].
- `strategy::CollapseStrategy=Hardest()`: Strategy for collapsing lineages.
  - `Hardest()` selects only the most frequent clone per lineage.
  - `Soft(cutoff)` keeps clones whose frequency within a lineage is at or above `cutoff`.
- `tie_breaker::AbstractTieBreaker=ByVdjCount()`: Tie-breaking policy when multiple clones
  share the maximum `clone_frequency` under `Hardest()`. Options: `ByVdjCount()`
  (default), `ByCdr3Count()`, `ByLexicographic()`, `BySequenceCount()`, `ByFirst()`,
  `ByMostNaive()`, `ByMostCommonVdjNt()`. Compose rules with `+` (e.g. `ByVdjCount() + ByMostNaive()`).
- `ByVdjCount()` requires a `vdj_nt` column (derived in `preprocess_data`
  when `v_sequence_start` and `j_sequence_end` are available).
- `ByMostNaive()` requires `v_identity` and `j_identity` columns.
- `ByMostCommonVdjNt()` matches igdiscover's clonotypes behavior: selects the row
  with the most common `vdj_nt` weighted by `count`. Requires `vdj_nt` column.
- `tie_atol::Real=0.0`: Absolute tolerance for considering clone frequencies equal
  when identifying ties.

# Returns
- `DataFrame`: Collapsed lineage data containing:
  - `clone_frequency`: Relative frequency of each clone within its lineage (0.0 to 1.0).
  - For `Hardest()` strategy:
    - `count`: Sum of counts across all lineage members
    - `nVDJ_nt`: Number of unique VDJ_nt sequences in the lineage

# Example
```julia
lineages = DataFrame(...)  # Your input data

# Default collapse (aggregates count, adds nVDJ_nt)
collapsed = collapse_lineages(lineages, Hardest())

# To match igdiscover's clonotypes output exactly:
collapsed = collapse_lineages(lineages, Hardest(); tie_breaker=ByMostCommonVdjNt())

# Soft collapse (keeps multiple clones per lineage)
collapsed = collapse_lineages(lineages, Soft(0.1))
```
"""
function collapse_lineages(df::DataFrame, strategy::CollapseStrategy=Hardest();
                           tie_breaker::AbstractTieBreaker=ByVdjCount(),
                           tie_atol::Real=0.0)
    with_frequency = clone_frequency_table(df)

    df_with_freq = leftjoin(df, with_frequency,
        on=[:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3],
        makeunique=true)

    if requires_vdj_nt(tie_breaker) && :vdj_nt ∉ propertynames(df_with_freq)
        throw(ArgumentError("vdj_nt column is required for the selected tie_breaker. " *
            "Provide vdj_nt or preprocess with v_sequence_start and j_sequence_end."))
    end

    if requires_vdj_count(tie_breaker)
        df_with_freq = add_vdj_count(df_with_freq)
    end

    aggregates = compute_aggregates(strategy, df_with_freq)

    collapsed = combine(groupby(df_with_freq, :lineage_id)) do group
        collapse_group(strategy, group, tie_breaker, tie_atol)
    end

    return apply_aggregates(collapsed, aggregates)
end

"""
    hardest_tie_summary(df::DataFrame; atol::Real=0.0)::DataFrame

Summarize lineages where multiple clones share the maximum clone frequency.
Returns a DataFrame with one row per lineage, including the count of tied clones
and the CDR3s involved in the tie. Set `atol` to a positive value to treat
frequencies within that tolerance as equal.
"""
function hardest_tie_summary(df::DataFrame; atol::Real=0.0)::DataFrame
    with_frequency = clone_frequency_table(df)
    return combine(groupby(with_frequency, :lineage_id)) do group
        max_freq = maximum(group.clone_frequency)
        if atol > 0
            tie_mask = abs.(group.clone_frequency .- max_freq) .<= atol
        else
            tie_mask = group.clone_frequency .== max_freq
        end
        tied_cdr3s = group.cdr3[tie_mask]
        (
            max_clone_frequency = max_freq,
            hardest_tie_count = length(tied_cdr3s),
            hardest_tied = length(tied_cdr3s) > 1,
            hardest_tied_cdr3s = [tied_cdr3s],
        )
    end
end