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
    n = nrow(df)
    
    # Build index arrays for computed columns - no copy of df
    cluster_vec = Vector{Int}(undef, n)
    min_dist_vec = Vector{Float32}(undef, n)
    cluster_size_vec = Vector{Int}(undef, n)
    cdr3_count_vec = Vector{Int}(undef, n)
    max_cdr3_count_vec = Vector{Int}(undef, n)
    cdr3_freq_vec = Vector{Float64}(undef, n)
    
    # Group by VJ+length and process each group
    groups = groupby(df, [:v_call_first, :j_call_first, :cdr3_length]; sort=false)
    group_indices = groupindices(groups)
    
    # Process groups in parallel, writing directly to output vectors
    ngroups = length(groups)
    
    # First pass: compute clusters per group
    Threads.@threads for gi in 1:ngroups
        group = groups[gi]
        rows = parentindices(group)[1]  # Get row indices in original df
        process_vj_group_inplace!(group, rows, cluster_vec, min_dist_vec, 
                                  cluster_size_vec, cdr3_count_vec, max_cdr3_count_vec, cdr3_freq_vec,
                                  distance_metric, clustering_method)
    end
    
    # Compute lineage IDs using groupby on original df + cluster
    # Build a temporary DataFrame with just the keys needed
    temp_keys = DataFrame(
        v_call_first = df.v_call_first,
        j_call_first = df.j_call_first, 
        cdr3_length = df.cdr3_length,
        cluster = cluster_vec
    )
    lineage_id_vec = groupindices(groupby(temp_keys, [:v_call_first, :j_call_first, :cdr3_length, :cluster]; sort=false))
    
    # Build result by combining original df with computed columns (no copy, uses views)
    result = hcat(df, DataFrame(
        cluster = cluster_vec,
        min_distance = min_dist_vec,
        cluster_size = cluster_size_vec,
        cdr3_count = cdr3_count_vec,
        max_cdr3_count = max_cdr3_count_vec,
        cdr3_frequency = cdr3_freq_vec,
        lineage_id = lineage_id_vec
    ); copycols=false)
    
    result
end

function process_vj_group_inplace!(group::SubDataFrame, rows::AbstractVector{Int},
                                   cluster_vec::Vector{Int}, min_dist_vec::Vector{Float32},
                                   cluster_size_vec::Vector{Int}, cdr3_count_vec::Vector{Int},
                                   max_cdr3_count_vec::Vector{Int}, cdr3_freq_vec::Vector{Float64},
                                   metric, clustering)
    cdr3_col = group.cdr3
    unique_cdr3 = unique(cdr3_col)
    n_unique = length(unique_cdr3)
    n_rows = length(rows)
    
    # Compute clusters and min distances for unique CDR3s
    local clusters::Vector{Int}
    local min_dists::Vector{Float32}
    
    if n_unique > 1
        seqs = LongDNA{4}.(unique_cdr3)
        dist = compute_pairwise_distance(metric, seqs)
        min_dists = compute_min_distances(dist)
        clusters = perform_clustering(clustering, :single, dist)
    else
        clusters = [1]
        min_dists = [0.0f0]
    end
    
    # Build lookup from CDR3 to cluster/min_dist
    cdr3_to_idx = Dict{eltype(unique_cdr3),Int}(unique_cdr3[i] => i for i in 1:n_unique)
    
    # Map to all rows and compute cluster stats
    cluster_sizes = Dict{Int,Int}()
    cluster_cdr3_counts = Dict{Tuple{Int,eltype(unique_cdr3)},Int}()
    
    @inbounds for (li, ri) in enumerate(rows)
        cdr3 = cdr3_col[li]
        idx = cdr3_to_idx[cdr3]
        c = clusters[idx]
        cluster_vec[ri] = c
        min_dist_vec[ri] = min_dists[idx]
        cluster_sizes[c] = get(cluster_sizes, c, 0) + 1
        cluster_cdr3_counts[(c, cdr3)] = get(cluster_cdr3_counts, (c, cdr3), 0) + 1
    end
    
    # Compute max cdr3 count per cluster
    cluster_max_cdr3 = Dict{Int,Int}()
    for ((c, _), cnt) in cluster_cdr3_counts
        cluster_max_cdr3[c] = max(get(cluster_max_cdr3, c, 0), cnt)
    end
    
    # Fill stats vectors
    @inbounds for (li, ri) in enumerate(rows)
        cdr3 = cdr3_col[li]
        c = cluster_vec[ri]
        cluster_size_vec[ri] = cluster_sizes[c]
        cnt = cluster_cdr3_counts[(c, cdr3)]
        cdr3_count_vec[ri] = cnt
        mx = cluster_max_cdr3[c]
        max_cdr3_count_vec[ri] = mx
        cdr3_freq_vec[ri] = cnt / mx
    end
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
    n = nrow(df)
    
    # Validate requirements
    if requires_vdj_nt(tie_breaker) && !hasproperty(df, :vdj_nt)
        throw(ArgumentError("vdj_nt column required for selected tie_breaker"))
    end
    
    # Compute clone frequencies into lookup table (no mutation of df)
    clone_freq, seq_count = compute_clone_frequency_lookup(df)
    
    # Compute vdj_count lookup if needed
    vdj_count_lookup = requires_vdj_count(tie_breaker) ? compute_vdj_count_lookup(df) : nothing
    
    # Compute aggregates before collapse (for Hardest strategy)  
    agg = compute_aggregates(strategy, df)
    
    # Select indices based on strategy
    selected_indices = collapse_by_strategy(df, strategy, tie_breaker, Float64(tie_atol), 
                                            clone_freq, seq_count, vdj_count_lookup)
    
    # Build result from selected rows
    build_collapse_result(df, selected_indices, strategy, clone_freq, seq_count, vdj_count_lookup, agg)
end

# Compute clone frequency as lookup table: (d_region, lineage_id, j_call, v_call, cdr3) -> (frequency, seq_count)
function compute_clone_frequency_lookup(df::DataFrame)
    clone_keys = [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3]
    clone_counts = combine(groupby(df, clone_keys; sort=false), nrow => :_cnt)
    lineage_totals = combine(groupby(clone_counts, :lineage_id; sort=false), :_cnt => sum => :_total)
    
    # Build lookup
    freq_lookup = Dict{Tuple{Any,Int,Any,Any,Any}, Float64}()
    count_lookup = Dict{Tuple{Any,Int,Any,Any,Any}, Int}()
    
    # Create lineage total lookup
    lin_totals = Dict{Int,Int}(row.lineage_id => row._total for row in eachrow(lineage_totals))
    
    for row in eachrow(clone_counts)
        key = (row.d_region, row.lineage_id, row.j_call_first, row.v_call_first, row.cdr3)
        freq_lookup[key] = row._cnt / lin_totals[row.lineage_id]
        count_lookup[key] = row._cnt
    end
    
    freq_lookup, count_lookup
end

function compute_vdj_count_lookup(df::DataFrame)
    has_count = hasproperty(df, :count)
    
    grouped = if !has_count
        combine(groupby(df, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3, :vdj_nt]; sort=false),
                nrow => :_vdj_cnt)
    else
        combine(groupby(df, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3, :vdj_nt]; sort=false),
                :count => sum => :_vdj_cnt)
    end
    
    # Max vdj_count per clone
    by_clone = combine(groupby(grouped, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3]; sort=false),
                       :_vdj_cnt => maximum => :_max_vdj)
    
    lookup = Dict{Tuple{Int,Any,Any,Any,Any}, Int}()
    for row in eachrow(by_clone)
        key = (row.lineage_id, row.d_region, row.j_call_first, row.v_call_first, row.cdr3)
        lookup[key] = row._max_vdj
    end
    lookup
end

function collapse_by_strategy(df::DataFrame, ::Hardest, tie_breaker::AbstractTieBreaker, atol::Float64,
                              clone_freq, seq_count, vdj_count_lookup)
    groups = groupby(df, :lineage_id; sort=false)
    selected = Vector{Int}(undef, length(groups))
    
    for (gi, group) in enumerate(groups)
        rows = parentindices(group)[1]
        selected[gi] = select_hardest_index(df, rows, tie_breaker, atol, clone_freq, seq_count, vdj_count_lookup)
    end
    selected
end

function collapse_by_strategy(df::DataFrame, strategy::Soft, ::AbstractTieBreaker, ::Float64,
                              clone_freq, seq_count, vdj_count_lookup)
    cutoff = strategy.cutoff
    indices = Int[]
    sizehint!(indices, nrow(df) ÷ 2)
    
    @inbounds for i in 1:nrow(df)
        key = (df.d_region[i], df.lineage_id[i], df.j_call_first[i], df.v_call_first[i], df.cdr3[i])
        clone_freq[key] >= cutoff && push!(indices, i)
    end
    indices
end

function select_hardest_index(df::DataFrame, rows::AbstractVector{Int}, 
                              tie_breaker::AbstractTieBreaker, atol::Float64,
                              clone_freq, seq_count, vdj_count_lookup)
    # Find max frequency among rows
    max_freq = 0.0
    @inbounds for ri in rows
        key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])
        f = clone_freq[key]
        f > max_freq && (max_freq = f)
    end
    
    # Find candidates (rows with max frequency)
    candidates = if atol > 0.0
        Int[ri for ri in rows if begin
            key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])
            abs(clone_freq[key] - max_freq) <= atol
        end]
    else
        Int[ri for ri in rows if begin
            key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])
            clone_freq[key] == max_freq
        end]
    end
    
    select_representative_index(tie_breaker, df, rows, candidates, seq_count, vdj_count_lookup)
end

# MostCommonVdjNtTieBreaker - matches igdiscover's behavior exactly
function select_representative_index(::MostCommonVdjNtTieBreaker, df::DataFrame, 
                                     rows::AbstractVector{Int}, ::Vector{Int}, ::Any, ::Any)
    n = length(rows)
    n <= 2 && return rows[1]
    
    # Process in original row order (rows are already in df order within the group)
    # For igdiscover parity, we need insertion order - rows is already sorted
    vdj_col = df.vdj_nt
    count_col = hasproperty(df, :count) ? df.count : nothing
    
    vdj_counts = Dict{String,Int}()
    vdj_first_row = Dict{String,Int}()
    
    for ri in rows
        @inbounds vdj = vdj_col[ri]
        ismissing(vdj) && continue
        @inbounds cnt = isnothing(count_col) ? 1 : (ismissing(count_col[ri]) ? 1 : count_col[ri])
        
        if haskey(vdj_counts, vdj)
            vdj_counts[vdj] += cnt
        else
            vdj_counts[vdj] = cnt
            vdj_first_row[vdj] = ri
        end
    end
    
    isempty(vdj_counts) && return rows[1]
    
    # Find VDJ with max count, using first occurrence for ties
    max_count = 0
    best_row = rows[1]
    
    for (vdj, cnt) in vdj_counts
        first_row = vdj_first_row[vdj]
        if cnt > max_count || (cnt == max_count && first_row < best_row)
            max_count = cnt
            best_row = first_row
        end
    end
    
    best_row
end

# TieBreaker - uses sorting criteria  
function select_representative_index(tb::TieBreaker, df::DataFrame, 
                                     rows::AbstractVector{Int}, candidates::Vector{Int}, seq_count, vdj_count_lookup)
    length(candidates) == 1 && return candidates[1]
    isempty(tb.criteria) && return candidates[1]
    
    # Validate columns exist (sequence_count and vdj_count come from lookups)
    for (col, _) in tb.criteria
        col in (:vdj_count, :sequence_count) && continue
        hasproperty(df, col) || throw(ArgumentError("Missing column: $col"))
    end
    
    best = candidates[1]
    @inbounds for i in 2:length(candidates)
        idx = candidates[i]
        is_better_row(df, idx, best, tb.criteria, seq_count, vdj_count_lookup) && (best = idx)
    end
    best
end

@inline function is_better_row(df, idx1::Int, idx2::Int, criteria::Vector{Pair{Symbol,Bool}}, seq_count, vdj_count_lookup)
    @inbounds for (col, desc) in criteria
        v1 = get_lookup_or_col(df, idx1, col, seq_count, vdj_count_lookup)
        v2 = get_lookup_or_col(df, idx2, col, seq_count, vdj_count_lookup)
        isequal(v1, v2) && continue
        ismissing(v1) && return false
        ismissing(v2) && return true
        return desc ? v1 > v2 : v1 < v2
    end
    false
end

@inline function get_lookup_or_col(df, ri, col, seq_count, vdj_count_lookup)
    if col == :vdj_count
        return get_vdj_count(df, ri, vdj_count_lookup)
    elseif col == :sequence_count
        key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])
        return get(seq_count, key, 0)
    else
        return df[ri, col]
    end
end

@inline function get_vdj_count(df, ri, lookup)
    isnothing(lookup) && return 0
    key = (df.lineage_id[ri], df.d_region[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])
    get(lookup, key, 0)
end

requires_vdj_nt(::MostCommonVdjNtTieBreaker) = true
requires_vdj_nt(tb::TieBreaker) = any(col == :vdj_count for (col, _) in tb.criteria)
requires_vdj_count(::MostCommonVdjNtTieBreaker) = false
requires_vdj_count(tb::TieBreaker) = any(col == :vdj_count for (col, _) in tb.criteria)

# Aggregation types
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
        cnt = isnothing(count_col) ? 1 : (ismissing(count_col[i]) ? 0 : count_col[i])
        count_sums[lid] = get(count_sums, lid, 0) + cnt
        
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

function build_collapse_result(df::DataFrame, selected_indices::Vector{Int}, strategy::Hardest,
                               clone_freq, seq_count, vdj_count_lookup, agg::LineageAggregates)
    # Slice original df (single allocation)
    result = df[selected_indices, :]
    
    # Add aggregated columns
    n = length(selected_indices)
    counts = Vector{Int}(undef, n)
    nvdj = Vector{Union{Missing,Int}}(undef, n)
    has_vdj = !isempty(agg.vdj_counts)
    
    @inbounds for i in 1:n
        lid = result.lineage_id[i]
        counts[i] = get(agg.count_sums, lid, 0)
        nvdj[i] = has_vdj ? get(agg.vdj_counts, lid, missing) : missing
    end
    
    result[!, :count] = counts
    result[!, :nVDJ_nt] = nvdj
    result
end

function build_collapse_result(df::DataFrame, selected_indices::Vector{Int}, ::Soft,
                               clone_freq, seq_count, vdj_count_lookup, ::NoAggregates)
    result = df[selected_indices, :]
    
    # Add clone_frequency and sequence_count columns
    n = length(selected_indices)
    freq_vec = Vector{Float64}(undef, n)
    count_vec = Vector{Int}(undef, n)
    
    @inbounds for i in 1:n
        ri = selected_indices[i]
        key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])
        freq_vec[i] = clone_freq[key]
        count_vec[i] = seq_count[key]
    end
    
    result[!, :clone_frequency] = freq_vec
    result[!, :sequence_count] = count_vec
    result
end

"""
    hardest_tie_summary(df::DataFrame; atol=0.0) -> DataFrame

Diagnostic function to identify lineages with tied maximum clone frequencies.
"""
function hardest_tie_summary(df::DataFrame; atol::Real=0.0)
    clone_freq, _ = compute_clone_frequency_lookup(df)
    
    combine(groupby(df, :lineage_id; sort=false)) do group
        rows = parentindices(group)[1]
        
        # Find max frequency
        mx = 0.0
        for ri in rows
            key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])
            f = clone_freq[key]
            f > mx && (mx = f)
        end
        
        # Find tied CDR3s
        tied = String[]
        for ri in rows
            key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])
            f = clone_freq[key]
            is_tied = atol > 0 ? abs(f - mx) <= atol : f == mx
            is_tied && !(df.cdr3[ri] in tied) && push!(tied, df.cdr3[ri])
        end
        
        (max_clone_frequency=mx, hardest_tie_count=length(tied), 
         hardest_tied=length(tied) > 1, hardest_tied_cdr3s=[tied])
    end
end
