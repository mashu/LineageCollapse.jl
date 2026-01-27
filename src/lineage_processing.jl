"""
Abstract type for all distance metrics used in sequence comparison.
"""
abstract type AbstractDistanceMetric end

"""
Abstract type for clustering methods.
"""
abstract type ClusteringMethod end

"""
Abstract type for lineage collapse strategies.
"""
abstract type CollapseStrategy end

"""
    HammingDistance <: AbstractDistanceMetric

Hamming distance metric - counts the number of mismatches between sequences.
"""
struct HammingDistance <: AbstractDistanceMetric end

"""
    NormalizedHammingDistance <: AbstractDistanceMetric

Normalized Hamming distance - mismatches divided by sequence length.
"""
struct NormalizedHammingDistance <: AbstractDistanceMetric end

"""
    LevenshteinDistance <: AbstractDistanceMetric

Levenshtein (edit) distance metric.
"""
struct LevenshteinDistance <: AbstractDistanceMetric end

"""
    Hardest <: CollapseStrategy

Collapse strategy that selects exactly one representative per lineage.
"""
struct Hardest <: CollapseStrategy end

"""
    Soft{T} <: CollapseStrategy

Collapse strategy that keeps all clones above a frequency threshold.
"""
struct Soft{T<:AbstractFloat} <: CollapseStrategy
    cutoff::T
    function Soft{T}(cutoff::T) where {T<:AbstractFloat}
        (0.0 <= cutoff <= 1.0) || throw(ArgumentError("cutoff must be between 0.0 and 1.0"))
        new{T}(cutoff)
    end
end

Soft(cutoff::T) where {T<:AbstractFloat} = Soft{T}(cutoff)
Soft(cutoff::Integer) = Soft(float(cutoff))

"""
    HierarchicalClustering{T} <: ClusteringMethod

Hierarchical clustering with a distance cutoff.
"""
struct HierarchicalClustering{T<:AbstractFloat} <: ClusteringMethod
    cutoff::T
end

HierarchicalClustering(cutoff::Integer) = HierarchicalClustering{Float32}(Float32(cutoff))

# Fixed key types - preprocessing normalizes all strings to String
const CloneKey = Tuple{String, Int, String, String, String}  # d_region, lineage_id, j_call_first, v_call_first, cdr3
const VdjKey = Tuple{Int, String, String, String, String}    # lineage_id, d_region, j_call_first, v_call_first, cdr3

"""
    compute_distance(metric::AbstractDistanceMetric, x::LongDNA{4}, y::LongDNA{4}) -> Float32

Compute distance between two DNA sequences using the specified metric.
"""
function compute_distance end

@inline compute_distance(::HammingDistance, x::LongDNA{4}, y::LongDNA{4}) = 
    Float32(mismatches(x, y))

@inline compute_distance(::NormalizedHammingDistance, x::LongDNA{4}, y::LongDNA{4}) = 
    Float32(mismatches(x, y)) / Float32(length(x))

@inline compute_distance(::LevenshteinDistance, x::LongDNA{4}, y::LongDNA{4}) = 
    Float32(evaluate(Levenshtein(), String(x), String(y)))

"""
    compute_pairwise_distance(metric, sequences) -> Matrix{Float32}

Compute pairwise distance matrix for DNA sequences.
"""
function compute_pairwise_distance(metric::M, sequences::AbstractVector{S}) where {
        M<:AbstractDistanceMetric, S<:LongSequence{DNAAlphabet{4}}}
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
    perform_clustering(method::HierarchicalClustering, linkage, dist_matrix) -> Vector{Int}

Perform hierarchical clustering and return cluster assignments.
"""
function perform_clustering(method::HierarchicalClustering, linkage::Symbol, 
                           dist_matrix::AbstractMatrix{<:AbstractFloat})
    cutree(hclust(dist_matrix, linkage=linkage, uplo=:U), h=method.cutoff)
end

# Dispatch on threshold type to select metric
select_metric(::Type{<:Integer}) = HammingDistance()
select_metric(::Type{<:AbstractFloat}) = NormalizedHammingDistance()

select_clustering(threshold::Integer) = HierarchicalClustering(Float32(threshold))
select_clustering(threshold::AbstractFloat) = HierarchicalClustering(Float32(threshold))

"""
    process_lineages(df::DataFrame, threshold; linkage=:single) -> DataFrame

Process sequences into lineages using CDR3 clustering.
Threshold type determines metric: Integer -> Hamming, Float -> Normalized.
"""
function process_lineages(df::DataFrame, threshold::T; linkage::Symbol=:single) where {T<:Union{Integer,AbstractFloat}}
    threshold >= 0 || throw(ArgumentError("threshold must be non-negative"))
    T <: AbstractFloat && !(0.0 <= threshold <= 1.0) && 
        throw(ArgumentError("float threshold must be between 0.0 and 1.0"))
    
    process_lineages_impl(df, select_metric(T), select_clustering(threshold), linkage)
end

"""
    process_lineages(df::DataFrame; distance_metric, clustering_method, linkage) -> DataFrame

Process sequences into lineages with explicit metric and clustering configuration.
"""
function process_lineages(df::DataFrame;
                          distance_metric::AbstractDistanceMetric=HammingDistance(),
                          clustering_method::ClusteringMethod=HierarchicalClustering(1.0f0),
                          linkage::Symbol=:single)
    process_lineages_impl(df, distance_metric, clustering_method, linkage)
end

function process_lineages_impl(df::DataFrame, metric::M, clustering::C, linkage::Symbol) where {
        M<:AbstractDistanceMetric, C<:ClusteringMethod}
    n = nrow(df)
    
    # Pre-allocate output vectors
    cluster_vec = Vector{Int}(undef, n)
    min_dist_vec = Vector{Float32}(undef, n)
    cluster_size_vec = Vector{Int}(undef, n)
    cdr3_count_vec = Vector{Int}(undef, n)
    max_cdr3_count_vec = Vector{Int}(undef, n)
    cdr3_freq_vec = Vector{Float64}(undef, n)
    
    groups = groupby(df, [:v_call_first, :j_call_first, :cdr3_length]; sort=false)
    ngroups = length(groups)
    
    Threads.@threads for gi in 1:ngroups
        group = groups[gi]
        rows = parentindices(group)[1]
        process_vj_group!(group, rows, cluster_vec, min_dist_vec, 
                         cluster_size_vec, cdr3_count_vec, max_cdr3_count_vec, cdr3_freq_vec,
                         metric, clustering, linkage)
    end
    
    # Compute lineage IDs
    temp_keys = DataFrame(
        v_call_first = df.v_call_first,
        j_call_first = df.j_call_first, 
        cdr3_length = df.cdr3_length,
        cluster = cluster_vec
    )
    lineage_id_vec = groupindices(groupby(temp_keys, [:v_call_first, :j_call_first, :cdr3_length, :cluster]; sort=false))
    
    hcat(df, DataFrame(
        cluster = cluster_vec,
        min_distance = min_dist_vec,
        cluster_size = cluster_size_vec,
        cdr3_count = cdr3_count_vec,
        max_cdr3_count = max_cdr3_count_vec,
        cdr3_frequency = cdr3_freq_vec,
        lineage_id = lineage_id_vec
    ); copycols=false)
end

function process_vj_group!(group::SubDataFrame, rows::AbstractVector{Int},
                          cluster_vec::Vector{Int}, min_dist_vec::Vector{Float32},
                          cluster_size_vec::Vector{Int}, cdr3_count_vec::Vector{Int},
                          max_cdr3_count_vec::Vector{Int}, cdr3_freq_vec::Vector{Float64},
                          metric::M, clustering::C, linkage::Symbol) where {M<:AbstractDistanceMetric, C<:ClusteringMethod}
    cdr3_col = group.cdr3
    unique_cdr3 = unique(cdr3_col)
    n_unique = length(unique_cdr3)
    
    local clusters::Vector{Int}
    local min_dists::Vector{Float32}
    
    if n_unique > 1
        seqs = LongDNA{4}.(unique_cdr3)
        dist = compute_pairwise_distance(metric, seqs)
        min_dists = compute_min_distances(dist)
        clusters = perform_clustering(clustering, linkage, dist)
    else
        clusters = Int[1]
        min_dists = Float32[0.0f0]
    end
    
    # Build lookup - use String keys since preprocessing normalizes types
    cdr3_to_idx = Dict{String,Int}(unique_cdr3[i] => i for i in 1:n_unique)
    cluster_sizes = Dict{Int,Int}()
    cluster_cdr3_counts = Dict{Tuple{Int,String},Int}()
    
    @inbounds for (li, ri) in enumerate(rows)
        cdr3 = cdr3_col[li]
        idx = cdr3_to_idx[cdr3]
        c = clusters[idx]
        cluster_vec[ri] = c
        min_dist_vec[ri] = min_dists[idx]
        cluster_sizes[c] = get(cluster_sizes, c, 0) + 1
        cluster_cdr3_counts[(c, cdr3)] = get(cluster_cdr3_counts, (c, cdr3), 0) + 1
    end
    
    cluster_max_cdr3 = Dict{Int,Int}()
    for ((c, _), cnt) in cluster_cdr3_counts
        cluster_max_cdr3[c] = max(get(cluster_max_cdr3, c, 0), cnt)
    end
    
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

# Lookup structures with fixed key types
struct CloneFrequencyLookup
    freq::Dict{CloneKey,Float64}
    count::Dict{CloneKey,Int}
end

struct VdjCountLookup
    lookup::Dict{VdjKey,Int}
end

struct NoVdjCountLookup end

"""
    collapse_lineages(df::DataFrame, strategy=Hardest(); tie_breaker, tie_atol) -> DataFrame

Collapse lineages to representative sequences.
"""
function collapse_lineages(df::DataFrame, strategy::CollapseStrategy=Hardest();
                           tie_breaker::AbstractTieBreaker=ByMostCommonVdjNt(),
                           tie_atol::Real=0.0)
    validate_tie_breaker(tie_breaker, df)
    
    clone_lookup = compute_clone_frequency_lookup(df)
    vdj_lookup = compute_vdj_count_lookup(tie_breaker, df)
    agg = compute_aggregates(strategy, df)
    
    selected = collapse_by_strategy(df, strategy, tie_breaker, Float64(tie_atol), clone_lookup, vdj_lookup)
    build_collapse_result(df, selected, strategy, clone_lookup, agg)
end

# Validation via dispatch
validate_tie_breaker(::MostCommonVdjNtTieBreaker, df::DataFrame) = 
    hasproperty(df, :vdj_nt) || throw(ArgumentError("vdj_nt column required"))

function validate_tie_breaker(tb::TieBreaker, df::DataFrame)
    any(col == :vdj_count for (col, _) in tb.criteria) && !hasproperty(df, :vdj_nt) && 
        throw(ArgumentError("vdj_nt column required for vdj_count"))
end

validate_tie_breaker(::AbstractTieBreaker, ::DataFrame) = nothing

# VDJ count lookup via dispatch
compute_vdj_count_lookup(::MostCommonVdjNtTieBreaker, ::DataFrame) = NoVdjCountLookup()

function compute_vdj_count_lookup(tb::TieBreaker, df::DataFrame)
    any(col == :vdj_count for (col, _) in tb.criteria) ? compute_vdj_count_lookup_impl(df) : NoVdjCountLookup()
end

function compute_clone_frequency_lookup(df::DataFrame)
    clone_keys = [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3]
    clone_counts = combine(groupby(df, clone_keys; sort=false), nrow => :_cnt)
    lineage_totals = combine(groupby(clone_counts, :lineage_id; sort=false), :_cnt => sum => :_total)
    
    freq = Dict{CloneKey,Float64}()
    cnt_dict = Dict{CloneKey,Int}()
    lin_totals = Dict{Int,Int}(row.lineage_id => row._total for row in eachrow(lineage_totals))
    
    for row in eachrow(clone_counts)
        key = (row.d_region, row.lineage_id, row.j_call_first, row.v_call_first, row.cdr3)::CloneKey
        freq[key] = row._cnt / lin_totals[row.lineage_id]
        cnt_dict[key] = row._cnt
    end
    
    CloneFrequencyLookup(freq, cnt_dict)
end

function compute_vdj_count_lookup_impl(df::DataFrame)
    has_count = hasproperty(df, :count)
    
    grouped = if !has_count
        combine(groupby(df, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3, :vdj_nt]; sort=false),
                nrow => :_vdj_cnt)
    else
        combine(groupby(df, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3, :vdj_nt]; sort=false),
                :count => sum => :_vdj_cnt)
    end
    
    by_clone = combine(groupby(grouped, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3]; sort=false),
                       :_vdj_cnt => maximum => :_max_vdj)
    
    lookup = Dict{VdjKey,Int}()
    for row in eachrow(by_clone)
        key = (row.lineage_id, row.d_region, row.j_call_first, row.v_call_first, row.cdr3)::VdjKey
        lookup[key] = row._max_vdj
    end
    
    VdjCountLookup(lookup)
end

# Collapse strategy dispatch
function collapse_by_strategy(df::DataFrame, ::Hardest, tie_breaker::AbstractTieBreaker, atol::Float64,
                              clone_lookup::CloneFrequencyLookup, vdj_lookup)
    groups = groupby(df, :lineage_id; sort=false)
    selected = Vector{Int}(undef, length(groups))
    
    for (gi, group) in enumerate(groups)
        rows = parentindices(group)[1]
        selected[gi] = select_hardest_index(df, rows, tie_breaker, atol, clone_lookup, vdj_lookup)
    end
    selected
end

function collapse_by_strategy(df::DataFrame, strategy::Soft, ::AbstractTieBreaker, ::Float64,
                              clone_lookup::CloneFrequencyLookup, vdj_lookup)
    cutoff = strategy.cutoff
    indices = Int[]
    sizehint!(indices, nrow(df) ÷ 2)
    
    @inbounds for i in 1:nrow(df)
        key = (df.d_region[i], df.lineage_id[i], df.j_call_first[i], df.v_call_first[i], df.cdr3[i])::CloneKey
        clone_lookup.freq[key] >= cutoff && push!(indices, i)
    end
    indices
end

function select_hardest_index(df::DataFrame, rows::AbstractVector{Int}, 
                              tie_breaker::AbstractTieBreaker, atol::Float64,
                              clone_lookup::CloneFrequencyLookup, vdj_lookup)
    max_freq = 0.0
    @inbounds for ri in rows
        key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])::CloneKey
        f = clone_lookup.freq[key]
        f > max_freq && (max_freq = f)
    end
    
    candidates = if atol > 0.0
        Int[ri for ri in rows if begin
            key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])::CloneKey
            abs(clone_lookup.freq[key] - max_freq) <= atol
        end]
    else
        Int[ri for ri in rows if begin
            key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])::CloneKey
            clone_lookup.freq[key] == max_freq
        end]
    end
    
    select_representative(tie_breaker, df, rows, candidates, clone_lookup, vdj_lookup)
end

# Representative selection via dispatch on tie_breaker type
function select_representative(::MostCommonVdjNtTieBreaker, df::DataFrame, 
                              rows::AbstractVector{Int}, ::Vector{Int}, 
                              ::CloneFrequencyLookup, ::NoVdjCountLookup)
    n = length(rows)
    n <= 2 && return rows[1]
    
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

function select_representative(tb::TieBreaker, df::DataFrame, 
                              ::AbstractVector{Int}, candidates::Vector{Int},
                              clone_lookup::CloneFrequencyLookup, vdj_lookup)
    length(candidates) == 1 && return candidates[1]
    isempty(tb.criteria) && return candidates[1]
    
    best = candidates[1]
    @inbounds for i in 2:length(candidates)
        idx = candidates[i]
        is_better(df, idx, best, tb.criteria, clone_lookup, vdj_lookup) && (best = idx)
    end
    best
end

# Comparison dispatch on lookup type
@inline function is_better(df, idx1::Int, idx2::Int, criteria::Vector{Pair{Symbol,Bool}},
                          clone_lookup::CloneFrequencyLookup, vdj_lookup::VdjCountLookup)
    @inbounds for (col, desc) in criteria
        v1 = get_value(df, idx1, col, clone_lookup, vdj_lookup)
        v2 = get_value(df, idx2, col, clone_lookup, vdj_lookup)
        isequal(v1, v2) && continue
        ismissing(v1) && return false
        ismissing(v2) && return true
        return desc ? v1 > v2 : v1 < v2
    end
    false
end

@inline function is_better(df, idx1::Int, idx2::Int, criteria::Vector{Pair{Symbol,Bool}},
                          clone_lookup::CloneFrequencyLookup, ::NoVdjCountLookup)
    @inbounds for (col, desc) in criteria
        v1 = get_value_simple(df, idx1, col, clone_lookup)
        v2 = get_value_simple(df, idx2, col, clone_lookup)
        isequal(v1, v2) && continue
        ismissing(v1) && return false
        ismissing(v2) && return true
        return desc ? v1 > v2 : v1 < v2
    end
    false
end

@inline function get_value(df, ri::Int, col::Symbol, clone_lookup::CloneFrequencyLookup, 
                          vdj_lookup::VdjCountLookup)
    if col === :sequence_count
        key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])::CloneKey
        return get(clone_lookup.count, key, 0)
    elseif col === :vdj_count
        key = (df.lineage_id[ri], df.d_region[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])::VdjKey
        return get(vdj_lookup.lookup, key, 0)
    else
        return df[ri, col]
    end
end

@inline function get_value_simple(df, ri::Int, col::Symbol, clone_lookup::CloneFrequencyLookup)
    if col === :sequence_count
        key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])::CloneKey
        return get(clone_lookup.count, key, 0)
    else
        return df[ri, col]
    end
end

# Requirement checks via dispatch
requires_vdj_nt(::MostCommonVdjNtTieBreaker) = true
requires_vdj_nt(tb::TieBreaker) = any(c === :vdj_count for (c, _) in tb.criteria)
requires_vdj_count(::MostCommonVdjNtTieBreaker) = false
requires_vdj_count(tb::TieBreaker) = any(c === :vdj_count for (c, _) in tb.criteria)

# Aggregation types
struct LineageAggregates
    count_sums::Dict{Int,Int}
    vdj_counts::Dict{Int,Int}
end

struct NoAggregates end

# Aggregation dispatch on strategy
function compute_aggregates(::Hardest, df::AbstractDataFrame)
    has_count = hasproperty(df, :count)
    has_vdj = hasproperty(df, :vdj_nt)
    
    count_sums = Dict{Int,Int}()
    vdj_sets = Dict{Int,Set{String}}()
    
    lineage_col = df.lineage_id
    
    @inbounds for i in 1:nrow(df)
        lid = lineage_col[i]
        cnt = has_count ? (ismissing(df.count[i]) ? 0 : df.count[i]) : 1
        count_sums[lid] = get(count_sums, lid, 0) + cnt
        
        if has_vdj
            vdj = df.vdj_nt[i]
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

# Result building dispatch on strategy
function build_collapse_result(df::DataFrame, selected::Vector{Int}, ::Hardest,
                               clone_lookup::CloneFrequencyLookup, agg::LineageAggregates)
    result = df[selected, :]
    n = length(selected)
    
    counts = Vector{Int}(undef, n)
    nvdj = Vector{Int}(undef, n)
    has_vdj = !isempty(agg.vdj_counts)
    
    @inbounds for i in 1:n
        lid = result.lineage_id[i]
        counts[i] = get(agg.count_sums, lid, 0)
        nvdj[i] = has_vdj ? get(agg.vdj_counts, lid, 0) : 0
    end
    
    result[!, :count] = counts
    result[!, :nVDJ_nt] = has_vdj ? nvdj : fill(missing, n)
    result
end

function build_collapse_result(df::DataFrame, selected::Vector{Int}, ::Soft,
                               clone_lookup::CloneFrequencyLookup, ::NoAggregates)
    result = df[selected, :]
    n = length(selected)
    
    freq_vec = Vector{Float64}(undef, n)
    count_vec = Vector{Int}(undef, n)
    
    @inbounds for i in 1:n
        ri = selected[i]
        key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])::CloneKey
        freq_vec[i] = clone_lookup.freq[key]
        count_vec[i] = clone_lookup.count[key]
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
    clone_lookup = compute_clone_frequency_lookup(df)
    
    combine(groupby(df, :lineage_id; sort=false)) do group
        rows = parentindices(group)[1]
        
        mx = 0.0
        for ri in rows
            key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])::CloneKey
            f = clone_lookup.freq[key]
            f > mx && (mx = f)
        end
        
        tied = String[]
        for ri in rows
            key = (df.d_region[ri], df.lineage_id[ri], df.j_call_first[ri], df.v_call_first[ri], df.cdr3[ri])::CloneKey
            f = clone_lookup.freq[key]
            is_tied = atol > 0 ? abs(f - mx) <= atol : f == mx
            is_tied && !(df.cdr3[ri] in tied) && push!(tied, df.cdr3[ri])
        end
        
        (max_clone_frequency=mx, hardest_tie_count=length(tied), 
         hardest_tied=length(tied) > 1, hardest_tied_cdr3s=[tied])
    end
end
