"""
Abstract type for distance metrics used in sequence comparison.
"""
abstract type DistanceMetric end

"""
Abstract type for normalized distance metrics (fraction of sequence length).
"""
abstract type NormalizedDistanceMetric end

"""
Abstract type for clustering methods.
"""
abstract type ClusteringMethod end

"""
Abstract type for lineage collapse strategies.
"""
abstract type CollapseStrategy end

"""
    HammingDistance <: DistanceMetric

Hamming distance metric - counts the number of mismatches between sequences.
"""
struct HammingDistance <: DistanceMetric end

"""
    NormalizedHammingDistance <: NormalizedDistanceMetric

Normalized Hamming distance - mismatches divided by sequence length.
"""
struct NormalizedHammingDistance <: NormalizedDistanceMetric end

"""
    LevenshteinDistance <: DistanceMetric

Levenshtein (edit) distance metric.
"""
struct LevenshteinDistance <: DistanceMetric end

"""
    Hardest <: CollapseStrategy

Collapse strategy that selects exactly one representative per lineage
(the clone with highest frequency within the lineage).
"""
struct Hardest <: CollapseStrategy end

"""
    Soft{T} <: CollapseStrategy

Collapse strategy that keeps all clones above a frequency threshold.

# Fields
- `cutoff::T`: Minimum clone frequency to retain (0.0 to 1.0)
"""
struct Soft{T<:AbstractFloat} <: CollapseStrategy
    cutoff::T
    function Soft{T}(cutoff::T) where {T<:AbstractFloat}
        (0.0 <= cutoff <= 1.0) || throw(ArgumentError("cutoff must be between 0.0 and 1.0"))
        new{T}(cutoff)
    end
end

Soft(cutoff::AbstractFloat) = Soft{typeof(cutoff)}(cutoff)
Soft(cutoff::Integer) = Soft(float(cutoff))

"""
    HierarchicalClustering <: ClusteringMethod

Hierarchical clustering with a distance cutoff.

# Fields
- `cutoff::Float32`: Distance threshold for cluster merging
"""
struct HierarchicalClustering <: ClusteringMethod
    cutoff::Float32
end

"""
    compute_distance(metric, x::LongDNA{4}, y::LongDNA{4}) -> Float32

Compute distance between two DNA sequences using the specified metric.
"""
function compute_distance end

@inline function compute_distance(::HammingDistance, x::LongDNA{4}, y::LongDNA{4})::Float32
    Float32(mismatches(x, y))
end

@inline function compute_distance(::NormalizedHammingDistance, x::LongDNA{4}, y::LongDNA{4})::Float32
    Float32(mismatches(x, y)) / Float32(length(x))
end

@inline function compute_distance(::LevenshteinDistance, x::LongDNA{4}, y::LongDNA{4})::Float32
    Float32(evaluate(Levenshtein(), String(x), String(y)))
end

"""
    compute_pairwise_distance(metric, sequences) -> Matrix{Float32}

Compute pairwise distance matrix for a collection of DNA sequences.
Returns upper triangular matrix.
"""
function compute_pairwise_distance(
    metric::M,
    sequences::AbstractVector{S}
)::Matrix{Float32} where {M<:Union{DistanceMetric,NormalizedDistanceMetric}, S<:LongSequence{DNAAlphabet{4}}}
    n = length(sequences)
    n > 0 || throw(ArgumentError("No sequences provided"))
    dist_matrix = zeros(Float32, n, n)

    if n < 50
        @inbounds for i in 1:n, j in (i+1):n
            dist_matrix[i, j] = compute_distance(metric, sequences[i], sequences[j])
        end
    else
        Threads.@threads for i in 1:n
            @inbounds for j in (i+1):n
                dist_matrix[i, j] = compute_distance(metric, sequences[i], sequences[j])
            end
        end
    end
    dist_matrix
end

"""
    perform_clustering(method::HierarchicalClustering, linkage::Symbol, dist_matrix) -> Vector{Int}

Perform hierarchical clustering and return cluster assignments.
"""
function perform_clustering(method::HierarchicalClustering, linkage::Symbol, dist_matrix::AbstractMatrix)::Vector{Int}
    cutree(hclust(dist_matrix, linkage=linkage, uplo=:U), h=method.cutoff)
end

@inline function cdr3_threshold_params(threshold::Integer)
    threshold >= 0 || throw(ArgumentError("threshold must be non-negative"))
    HammingDistance(), HierarchicalClustering(Float32(threshold))
end

@inline function cdr3_threshold_params(threshold::AbstractFloat)
    (0.0 <= threshold <= 1.0) || throw(ArgumentError("threshold must be between 0.0 and 1.0"))
    NormalizedHammingDistance(), HierarchicalClustering(Float32(threshold))
end

"""
    process_lineages(df::DataFrame, threshold; linkage=:single) -> DataFrame

Process sequences into lineages using CDR3 clustering.

# Arguments
- `df::DataFrame`: Input data with columns `v_call_first`, `j_call_first`, `cdr3`, `cdr3_length`
- `threshold`: If Integer, absolute mismatches. If Float, fraction of CDR3 length.
- `linkage::Symbol=:single`: Hierarchical clustering linkage method

# Returns
DataFrame with added columns: `lineage_id`, `cluster`, `min_distance`, `cluster_size`,
`cdr3_count`, `max_cdr3_count`, `cdr3_frequency`
"""
function process_lineages(df::DataFrame, threshold::Union{Integer,AbstractFloat}; linkage::Symbol=:single)::DataFrame
    metric, clustering = cdr3_threshold_params(threshold)
    process_lineages(df; distance_metric=metric, clustering_method=clustering, linkage=linkage)
end

"""
    process_lineages(df::DataFrame; distance_metric, clustering_method, linkage) -> DataFrame

Process sequences into lineages with explicit metric and clustering configuration.
"""
function process_lineages(df::DataFrame;
                          distance_metric::Union{DistanceMetric,NormalizedDistanceMetric}=HammingDistance(),
                          clustering_method::ClusteringMethod=HierarchicalClustering(1.0f0),
                          linkage::Symbol=:single)::DataFrame
    # Work on a copy to avoid mutating input
    work_df = copy(df)
    n = nrow(work_df)
    work_df[!, :_row_idx] = 1:n
    
    groups = groupby(work_df, [:v_call_first, :j_call_first, :cdr3_length]; sort=false)
    ngroups = length(groups)
    processed = Vector{DataFrame}(undef, ngroups)

    Threads.@threads for gi in 1:ngroups
        processed[gi] = process_vj_group(groups[gi], distance_metric, clustering_method, linkage)
    end

    result = reduce(vcat, processed)
    result[!, :lineage_id] = groupindices(groupby(result, [:v_call_first, :j_call_first, :cdr3_length, :cluster]; sort=false))
    
    sort!(result, :_row_idx)
    select!(result, Not(:_row_idx))
end

function process_vj_group(group::SubDataFrame, metric, clustering, linkage)
    unique_cdr3 = unique(group.cdr3)
    n_unique = length(unique_cdr3)
    group_df = DataFrame(group; copycols=false)

    if n_unique > 1
        seqs = LongDNA{4}.(unique_cdr3)
        dist = compute_pairwise_distance(metric, seqs)
        min_dists = compute_min_distances(dist)
        clusters = perform_clustering(clustering, linkage, dist)
        
        cdr3_to_cluster = Dict{eltype(unique_cdr3),Int}(unique_cdr3[i] => clusters[i] for i in 1:n_unique)
        cdr3_to_mindist = Dict{eltype(unique_cdr3),Float32}(unique_cdr3[i] => min_dists[i] for i in 1:n_unique)
        
        group_df[!, :cluster] = [cdr3_to_cluster[c] for c in group_df.cdr3]
        group_df[!, :min_distance] = [cdr3_to_mindist[c] for c in group_df.cdr3]
    else
        group_df[!, :cluster] = fill(1, nrow(group_df))
        group_df[!, :min_distance] = fill(0.0f0, nrow(group_df))
    end

    add_cluster_stats!(group_df)
end

function compute_min_distances(dist::Matrix{Float32})
    n = size(dist, 1)
    mins = fill(typemax(Float32), n)
    @inbounds for i in 1:n, j in (i+1):n
        d = dist[i, j]
        if d > 0.0f0
            d < mins[i] && (mins[i] = d)
            d < mins[j] && (mins[j] = d)
        end
    end
    @inbounds for i in 1:n
        mins[i] == typemax(Float32) && (mins[i] = 0.0f0)
    end
    mins
end

function add_cluster_stats!(df::DataFrame)
    n = nrow(df)
    cluster_col = df.cluster
    cdr3_col = df.cdr3
    
    cluster_sizes = Dict{Int,Int}()
    cluster_cdr3_counts = Dict{Tuple{Int,eltype(cdr3_col)},Int}()
    
    @inbounds for i in 1:n
        c = cluster_col[i]
        cdr3 = cdr3_col[i]
        cluster_sizes[c] = get(cluster_sizes, c, 0) + 1
        cluster_cdr3_counts[(c, cdr3)] = get(cluster_cdr3_counts, (c, cdr3), 0) + 1
    end
    
    cluster_max_cdr3 = Dict{Int,Int}()
    for ((c, _), cnt) in cluster_cdr3_counts
        cluster_max_cdr3[c] = max(get(cluster_max_cdr3, c, 0), cnt)
    end
    
    cluster_size_vec = Vector{Int}(undef, n)
    cdr3_count_vec = Vector{Int}(undef, n)
    max_cdr3_count_vec = Vector{Int}(undef, n)
    cdr3_freq_vec = Vector{Float64}(undef, n)
    
    @inbounds for i in 1:n
        c = cluster_col[i]
        cdr3 = cdr3_col[i]
        cluster_size_vec[i] = cluster_sizes[c]
        cnt = cluster_cdr3_counts[(c, cdr3)]
        cdr3_count_vec[i] = cnt
        mx = cluster_max_cdr3[c]
        max_cdr3_count_vec[i] = mx
        cdr3_freq_vec[i] = cnt / mx
    end
    
    df[!, :cluster_size] = cluster_size_vec
    df[!, :cdr3_count] = cdr3_count_vec
    df[!, :max_cdr3_count] = max_cdr3_count_vec
    df[!, :cdr3_frequency] = cdr3_freq_vec
    df
end

# Clone frequency computation - computes frequency of each clone within its lineage
function compute_clone_frequencies!(df::DataFrame)
    # Group by clone identity and count sequences
    clone_keys = [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3]
    clone_counts = combine(groupby(df, clone_keys; sort=false), nrow => :sequence_count)
    
    # Compute lineage totals
    lineage_totals = combine(groupby(clone_counts, :lineage_id; sort=false), :sequence_count => sum => :_lineage_total)
    
    # Join back and compute frequency
    clone_counts = leftjoin(clone_counts, lineage_totals, on=:lineage_id)
    clone_counts[!, :clone_frequency] = clone_counts.sequence_count ./ clone_counts._lineage_total
    select!(clone_counts, Not(:_lineage_total))
    
    # Join to original df
    leftjoin!(df, clone_counts, on=clone_keys, matchmissing=:equal)
    df
end

"""
    collapse_lineages(df::DataFrame, strategy=Hardest(); tie_breaker, tie_atol) -> DataFrame

Collapse lineages to representative sequences.

# Arguments
- `df::DataFrame`: Input data with `lineage_id` column from `process_lineages`
- `strategy::CollapseStrategy=Hardest()`: `Hardest()` for one per lineage, `Soft(cutoff)` for frequency threshold
- `tie_breaker::AbstractTieBreaker=ByMostCommonVdjNt()`: Tie-breaking strategy (matches igdiscover)
- `tie_atol::Real=0.0`: Tolerance for frequency comparison

# Returns
DataFrame with collapsed lineages. For `Hardest()`, includes aggregated `count` and `nVDJ_nt`.
"""
function collapse_lineages(df::DataFrame, strategy::CollapseStrategy=Hardest();
                           tie_breaker::AbstractTieBreaker=ByMostCommonVdjNt(),
                           tie_atol::Real=0.0)
    # Work on a copy to avoid mutating input
    work_df = copy(df)
    n = nrow(work_df)
    
    # Add row index for order preservation and final slicing
    work_df[!, :_row_idx] = 1:n
    
    # Compute clone frequencies (adds clone_frequency column)
    compute_clone_frequencies!(work_df)
    
    # Validate requirements
    if requires_vdj_nt(tie_breaker) && !hasproperty(work_df, :vdj_nt)
        throw(ArgumentError("vdj_nt column required for selected tie_breaker"))
    end
    
    # Add vdj_count if needed
    requires_vdj_count(tie_breaker) && add_vdj_count!(work_df)
    
    # Compute aggregates before collapse (for Hardest strategy)
    agg = compute_aggregates(strategy, work_df)
    
    # Collapse each lineage - returns vector of selected row indices
    selected_indices = collapse_by_strategy(work_df, strategy, tie_breaker, Float64(tie_atol))
    
    # Single slice at the end
    result = work_df[selected_indices, :]
    select!(result, Not(:_row_idx))
    
    # For Hardest, remove internal columns (replaced by aggregates)
    # For Soft, keep clone_frequency in output
    if strategy isa Hardest
        select!(result, Not([:clone_frequency, :sequence_count]))
    end
    
    # Apply aggregates
    apply_aggregates!(result, agg)
end

# Dispatch on strategy for collapse logic
function collapse_by_strategy(df::DataFrame, ::Hardest, tie_breaker::AbstractTieBreaker, atol::Float64)
    groups = groupby(df, :lineage_id; sort=false)
    selected = Vector{Int}(undef, length(groups))
    
    for (gi, group) in enumerate(groups)
        selected[gi] = select_hardest_index(group, tie_breaker, atol)
    end
    selected
end

function collapse_by_strategy(df::DataFrame, strategy::Soft, ::AbstractTieBreaker, ::Float64)
    cutoff = strategy.cutoff
    freq = df.clone_frequency
    indices = Int[]
    sizehint!(indices, nrow(df) ÷ 2)
    
    @inbounds for i in 1:nrow(df)
        freq[i] >= cutoff && push!(indices, i)
    end
    indices
end

# Select the best row index from a lineage group (Hardest strategy)
function select_hardest_index(group::SubDataFrame, tie_breaker::AbstractTieBreaker, atol::Float64)
    freq = group.clone_frequency
    max_freq = maximum(freq)
    n = nrow(group)
    
    # Find candidates (rows with max frequency)
    candidates = if atol > 0.0
        Int[i for i in 1:n if @inbounds abs(freq[i] - max_freq) <= atol]
    else
        Int[i for i in 1:n if @inbounds freq[i] == max_freq]
    end
    
    # Use tie-breaker to select one, then convert to global index
    local_idx = select_representative_index(tie_breaker, group, candidates)
    row_idx_col = group[!, :_row_idx]
    row_idx_col[local_idx]
end

# MostCommonVdjNtTieBreaker - matches igdiscover's behavior exactly
function select_representative_index(::MostCommonVdjNtTieBreaker, group::SubDataFrame, ::Vector{Int})
    n = nrow(group)
    n <= 2 && return 1
    
    # Sort by original order for consistent tie-breaking
    row_idx = group[!, :_row_idx]
    sorted_order = sortperm(row_idx)
    
    vdj_col = group.vdj_nt
    count_col = hasproperty(group, :count) ? group.count : nothing
    
    # Track VDJ_nt counts and first occurrence (in sorted order)
    vdj_counts = Dict{String,Int}()
    vdj_first = Dict{String,Int}()
    
    for (order_pos, i) in enumerate(sorted_order)
        @inbounds vdj = vdj_col[i]
        ismissing(vdj) && continue
        @inbounds cnt = isnothing(count_col) ? 1 : (ismissing(count_col[i]) ? 1 : count_col[i])
        
        if haskey(vdj_counts, vdj)
            vdj_counts[vdj] += cnt
        else
            vdj_counts[vdj] = cnt
            vdj_first[vdj] = order_pos
        end
    end
    
    isempty(vdj_counts) && return 1
    
    # Find VDJ with max count, using first occurrence for ties
    max_count = 0
    best_vdj = ""
    best_order = typemax(Int)
    
    for (vdj, cnt) in vdj_counts
        order = vdj_first[vdj]
        if cnt > max_count || (cnt == max_count && order < best_order)
            max_count = cnt
            best_vdj = vdj
            best_order = order
        end
    end
    
    # Find first row with this VDJ_nt (in sorted order)
    for i in sorted_order
        @inbounds if !ismissing(vdj_col[i]) && vdj_col[i] == best_vdj
            return i
        end
    end
    1
end

# TieBreaker - uses sorting criteria
function select_representative_index(tb::TieBreaker, group::SubDataFrame, candidates::Vector{Int})
    length(candidates) == 1 && return candidates[1]
    isempty(tb.criteria) && return candidates[1]
    
    # Validate columns exist
    for (col, _) in tb.criteria
        hasproperty(group, col) || throw(ArgumentError("Missing column: $col"))
    end
    
    best = candidates[1]
    @inbounds for i in 2:length(candidates)
        idx = candidates[i]
        is_better(group, idx, best, tb.criteria) && (best = idx)
    end
    best
end

@inline function is_better(df, idx1::Int, idx2::Int, criteria::Vector{Pair{Symbol,Bool}})
    @inbounds for (col, desc) in criteria
        v1, v2 = df[idx1, col], df[idx2, col]
        isequal(v1, v2) && continue
        ismissing(v1) && return false
        ismissing(v2) && return true
        return desc ? v1 > v2 : v1 < v2
    end
    false
end

requires_vdj_nt(::MostCommonVdjNtTieBreaker) = true
requires_vdj_nt(tb::TieBreaker) = any(col == :vdj_count for (col, _) in tb.criteria)
requires_vdj_count(::MostCommonVdjNtTieBreaker) = false
requires_vdj_count(tb::TieBreaker) = any(col == :vdj_count for (col, _) in tb.criteria)

# Add vdj_count column
function add_vdj_count!(df::DataFrame)
    has_count = hasproperty(df, :count)
    
    grouped = if !has_count
        combine(groupby(df, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3, :vdj_nt]; sort=false),
                nrow => :vdj_count)
    else
        combine(groupby(df, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3, :vdj_nt]; sort=false),
                :count => sum => :vdj_count)
    end

    by_clone = combine(groupby(grouped, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3]; sort=false),
                       :vdj_count => maximum => :vdj_count)

    leftjoin!(df, by_clone, on=[:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3], matchmissing=:equal)
    df
end

# Aggregation types for type-stable dispatch
struct LineageAggregates
    count_sums::Dict{Int,Int}
    vdj_counts::Dict{Int,Int}
end

struct NoAggregates end

function compute_aggregates(::Hardest, df::AbstractDataFrame)
    has_count = hasproperty(df, :count)
    has_vdj = hasproperty(df, :vdj_nt)
    
    count_sums = Dict{Int,Int}()
    vdj_sets = Dict{Int,Set{String}}()
    
    lineage_col = df.lineage_id
    count_col = has_count ? df.count : nothing
    vdj_col = has_vdj ? df.vdj_nt : nothing
    
    @inbounds for i in 1:nrow(df)
        lid = lineage_col[i]
        
        # Sum counts
        cnt = isnothing(count_col) ? 1 : (ismissing(count_col[i]) ? 0 : count_col[i])
        count_sums[lid] = get(count_sums, lid, 0) + cnt
        
        # Collect unique vdj_nt
        if !isnothing(vdj_col)
            vdj = vdj_col[i]
            if !ismissing(vdj)
                if !haskey(vdj_sets, lid)
                    vdj_sets[lid] = Set{String}()
                end
                push!(vdj_sets[lid], vdj)
            end
        end
    end
    
    vdj_counts = Dict{Int,Int}(lid => length(s) for (lid, s) in vdj_sets)
    LineageAggregates(count_sums, vdj_counts)
end

compute_aggregates(::Soft, ::AbstractDataFrame) = NoAggregates()

function apply_aggregates!(df::DataFrame, agg::LineageAggregates)
    n = nrow(df)
    counts = Vector{Int}(undef, n)
    nvdj = Vector{Union{Missing,Int}}(undef, n)
    
    has_vdj = !isempty(agg.vdj_counts)
    
    @inbounds for i in 1:n
        lid = df.lineage_id[i]
        counts[i] = get(agg.count_sums, lid, 0)
        nvdj[i] = has_vdj ? get(agg.vdj_counts, lid, missing) : missing
    end
    
    df[!, :count] = counts
    df[!, :nVDJ_nt] = nvdj
    df
end

apply_aggregates!(df::DataFrame, ::NoAggregates) = df

"""
    hardest_tie_summary(df::DataFrame; atol=0.0) -> DataFrame

Diagnostic function to identify lineages with tied maximum clone frequencies.
"""
function hardest_tie_summary(df::DataFrame; atol::Real=0.0)
    df_copy = copy(df)
    compute_clone_frequencies!(df_copy)
    
    combine(groupby(df_copy, :lineage_id; sort=false)) do group
        mx = maximum(group.clone_frequency)
        tied = atol > 0 ? 
            [group.cdr3[i] for i in 1:nrow(group) if abs(group.clone_frequency[i] - mx) <= atol] :
            [group.cdr3[i] for i in 1:nrow(group) if group.clone_frequency[i] == mx]
        (max_clone_frequency=mx, hardest_tie_count=length(tied), 
         hardest_tied=length(tied) > 1, hardest_tied_cdr3s=[tied])
    end
end
