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

Hamming distance metric — counts the number of mismatches between equal-length sequences.
"""
struct HammingDistance <: AbstractDistanceMetric end

"""
    NormalizedHammingDistance <: AbstractDistanceMetric

Normalized Hamming distance — mismatches divided by sequence length.
Requires equal-length sequences.
"""
struct NormalizedHammingDistance <: AbstractDistanceMetric end

"""
    LevenshteinDistance <: AbstractDistanceMetric

Levenshtein (edit) distance metric.

`process_lineages` groups by exact CDR3 length before clustering, so sequences of
different lengths are never compared. Levenshtein can still differ from Hamming for
equal-length strings (shifts/indels), but it cannot join CDR3s that differ in length.
"""
struct LevenshteinDistance <: AbstractDistanceMetric end

"""
    Hardest <: CollapseStrategy

Collapse strategy that selects exactly one representative per lineage from the
clone(s) with **maximum clone frequency**, then applies `tie_breaker` to those candidates.

Lineage-wide policies (`ByLineageMostCommonVdjNt`, `ByLowestVShm`) ignore the
frequency filter and select from every member of the lineage.
"""
struct Hardest <: CollapseStrategy end

"""
    IgDiscover <: CollapseStrategy

One representative per lineage using IgDiscover's `clonotypes.representative()`:

1. If the lineage has ≤ 2 rows, take the first row (input-order dependent).
2. Otherwise, take the first row whose `vdj_nt` has the largest lineage-wide
   summed abundance.

Clustering should use Hamming distance, threshold 1, and single linkage to match
IgDiscover's CDR3 grouping. Single linkage means sequences can share a lineage
through a chain of pairwise distances ≤ threshold even if some pairs exceed it.
"""
struct IgDiscover <: CollapseStrategy end

"""
    Soft{T} <: CollapseStrategy

Collapse strategy that keeps all clones above a frequency threshold.

# Fields
- `cutoff::T`: Minimum clone frequency to retain (0.0 to 1.0)

# Example
```julia
Soft(0.2)  # Keep clones with frequency ≥ 20%
```
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

# Fields
- `cutoff::T`: Distance threshold for cluster merging

# Example
```julia
HierarchicalClustering(1.0f0)  # Merge clusters within distance 1.0
```
"""
struct HierarchicalClustering{T<:AbstractFloat} <: ClusteringMethod
    cutoff::T
end

HierarchicalClustering(cutoff::Integer) = HierarchicalClustering{Float32}(Float32(cutoff))

# Fixed key types - preprocessing normalizes all strings to String
const CloneKey = Tuple{String, Int, String, String, String}  # d_region, lineage_id, j_call_first, v_call_first, cdr3
const VdjKey = Tuple{Int, String, String, String, String}    # lineage_id, d_region, j_call_first, v_call_first, cdr3
const DENSE_CLUSTER_WARN_N = 10_000

"""
    compute_distance(metric::AbstractDistanceMetric, x::LongDNA{4}, y::LongDNA{4}) -> Float32

Compute distance between two DNA sequences using the specified metric.
"""
function compute_distance end

function require_equal_length(x, y, metric)
    length(x) == length(y) || throw(ArgumentError(
        "$(typeof(metric)) requires equal-length sequences, got $(length(x)) and $(length(y))"
    ))
    return nothing
end

@inline function compute_distance(metric::HammingDistance, x::LongDNA{4}, y::LongDNA{4})
    require_equal_length(x, y, metric)
    Float32(mismatches(x, y))
end

@inline function compute_distance(metric::NormalizedHammingDistance, x::LongDNA{4}, y::LongDNA{4})
    require_equal_length(x, y, metric)
    Float32(mismatches(x, y)) / Float32(length(x))
end

@inline compute_distance(::LevenshteinDistance, x::LongDNA{4}, y::LongDNA{4}) =
    Float32(evaluate(Levenshtein(), String(x), String(y)))

"""
    compute_pairwise_distance(metric, sequences) -> Matrix{Float32}

Compute pairwise distance matrix for DNA sequences. Only the strict upper triangle
is filled; the lower triangle is left as zero for `hclust(..., uplo=:U)`.
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
`dist_matrix` may be a full symmetric matrix or an upper-triangle matrix used with `uplo=:U`.
"""
function perform_clustering(method::HierarchicalClustering, linkage::Symbol,
                           dist_matrix::AbstractMatrix{<:AbstractFloat})
    cutree(hclust(dist_matrix, linkage=linkage, uplo=:U), h=method.cutoff)
end

"""
    cluster_single_linkage(metric, sequences, cutoff) -> (clusters, min_dists)

Connected-component (single-linkage) clustering with O(n) extra memory.
Equivalent to IgDiscover's `hamming_single_linkage` when `metric` is Hamming and
`cutoff` is the mismatch threshold.
"""
function cluster_single_linkage(metric::M, sequences, cutoff::Real) where {M<:AbstractDistanceMetric}
    n = length(sequences)
    parent = collect(1:n)
    cutoff32 = Float32(cutoff)

    function find_root(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    function unite(i, j)
        ri = find_root(i)
        rj = find_root(j)
        ri == rj && return
        parent[rj] = ri
        return nothing
    end

    mins = fill(typemax(Float32), n)
    @inbounds for i in 1:n, j in (i+1):n
        d = compute_distance(metric, sequences[i], sequences[j])
        if d > 0.0f0
            d < mins[i] && (mins[i] = d)
            d < mins[j] && (mins[j] = d)
        end
        d <= cutoff32 && unite(i, j)
    end
    @inbounds for i in 1:n
        mins[i] == typemax(Float32) && (mins[i] = 0.0f0)
    end

    clusters = Vector{Int}(undef, n)
    idmap = Dict{Int,Int}()
    next_id = 1
    @inbounds for i in 1:n
        r = find_root(i)
        cid = get(idmap, r, 0)
        if cid == 0
            cid = next_id
            idmap[r] = next_id
            next_id += 1
        end
        clusters[i] = cid
    end
    return clusters, mins
end

# Dispatch on threshold type to select metric
select_metric(::Type{<:Integer}) = HammingDistance()
select_metric(::Type{<:AbstractFloat}) = NormalizedHammingDistance()

select_clustering(threshold::Integer) = HierarchicalClustering(Float32(threshold))
select_clustering(threshold::AbstractFloat) = HierarchicalClustering(Float32(threshold))

uses_frequency_candidates(::MostCommonVdjNtTieBreaker) = true
uses_frequency_candidates(::LineageMostCommonVdjNtTieBreaker) = false
uses_frequency_candidates(::LowestVShmTieBreaker) = false
uses_frequency_candidates(::AbstractTieBreaker) = true

default_tie_breaker(::Hardest) = ByMostCommonVdjNt()
default_tie_breaker(::Soft) = ByMostCommonVdjNt()
default_tie_breaker(::IgDiscover) = ByLineageMostCommonVdjNt(; small_lineage_first_row=true)

function validate_abundance(abundance::Symbol)
    (abundance === :count || abundance === :rows) ||
        throw(ArgumentError("abundance must be :count (use the count column when present) or :rows (each row is 1)"))
    return abundance
end

@inline function row_abundance(df::AbstractDataFrame, i::Int, abundance::Symbol)
    abundance === :rows && return 1
    hasproperty(df, :count) || return 1
    c = df.count[i]
    return ismissing(c) ? 1 : Int(c)
end

"""
    process_lineages(df::DataFrame, threshold; linkage=:single) -> DataFrame

Process sequences into lineages using CDR3 clustering.

Sequences are grouped by `v_call_first`, `j_call_first`, and `cdr3_length`. Within each
group, CDR3s are clustered. The default **single** linkage puts A and C in the same
lineage whenever there is a chain A—B—C of pairwise distances ≤ threshold, even if
A and C themselves exceed the threshold. Use `:complete` if every pair must be within
the threshold.

# Arguments
- `df::DataFrame`: Preprocessed data with columns `v_call_first`, `j_call_first`, `cdr3`, `cdr3_length`, `d_region`
- `threshold`: Clustering threshold. Integer for absolute mismatches, Float (0.0-1.0) for fraction of CDR3 length.
- `linkage::Symbol=:single`: Hierarchical clustering linkage (`:single`, `:complete`, `:average`)

# Returns
DataFrame with added columns:
- `lineage_id::Int`: Unique lineage identifier
- `cluster::Int`: Cluster assignment within V/J group
- `cluster_size::Int`: Number of sequences in cluster
- `min_distance::Float32`: Minimum distance to other CDR3s in cluster
- `cdr3_count::Int`: Count of this CDR3 in cluster
- `max_cdr3_count::Int`: Maximum CDR3 count in cluster
- `relative_to_max_cdr3::Float64`: `cdr3_count / max_cdr3_count` (not a conventional frequency)
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

# Keyword Arguments
- `distance_metric::AbstractDistanceMetric=HammingDistance()`: Distance metric for CDR3 comparison
- `clustering_method::ClusteringMethod=HierarchicalClustering(1.0f0)`: Clustering method and cutoff
- `linkage::Symbol=:single`: Hierarchical clustering linkage

See `process_lineages(df, threshold)` for return value documentation.
"""
function process_lineages(df::DataFrame;
                          distance_metric::AbstractDistanceMetric=HammingDistance(),
                          clustering_method::ClusteringMethod=HierarchicalClustering(1.0f0),
                          linkage::Symbol=:single)
    process_lineages_impl(df, distance_metric, clustering_method, linkage)
end

function empty_processed_lineages(df::DataFrame)
    hcat(df, DataFrame(
        cluster = Int[],
        min_distance = Float32[],
        cluster_size = Int[],
        cdr3_count = Int[],
        max_cdr3_count = Int[],
        relative_to_max_cdr3 = Float64[],
        lineage_id = Int[],
    ); copycols=false)
end

function process_lineages_impl(df::DataFrame, metric::M, clustering::C, linkage::Symbol) where {
        M<:AbstractDistanceMetric, C<:ClusteringMethod}
    n = nrow(df)
    n == 0 && return empty_processed_lineages(df)

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
        relative_to_max_cdr3 = cdr3_freq_vec,
        lineage_id = lineage_id_vec
    ); copycols=false)
end

function cluster_unique_sequences(metric, seqs, clustering::HierarchicalClustering, ::Val{:single})
    cluster_single_linkage(metric, seqs, clustering.cutoff)
end

function cluster_unique_sequences(metric, seqs, clustering::HierarchicalClustering, ::Val{L}) where L
    n = length(seqs)
    n >= DENSE_CLUSTER_WARN_N && @warn "Dense hierarchical clustering of $n unique CDR3s allocates an n×n Float32 distance matrix"
    dist = compute_pairwise_distance(metric, seqs)
    min_dists = compute_min_distances(dist)
    clusters = perform_clustering(clustering, L, dist)
    return clusters, min_dists
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
        clusters, min_dists = cluster_unique_sequences(metric, seqs, clustering, Val(linkage))
    else
        clusters = Int[1]
        min_dists = Float32[0.0f0]
    end

    n_clusters = maximum(clusters)
    cdr3_to_idx = Dict{String,Int}(unique_cdr3[i] => i for i in 1:n_unique)

    # Use Vectors for cluster-indexed data (clusters are 1:K)
    cluster_sizes = zeros(Int, n_clusters)
    cluster_max_cdr3 = zeros(Int, n_clusters)
    cluster_cdr3_counts = Dict{Tuple{Int,String},Int}()

    @inbounds for (li, ri) in enumerate(rows)
        cdr3 = cdr3_col[li]
        idx = cdr3_to_idx[cdr3]
        c = clusters[idx]
        cluster_vec[ri] = c
        min_dist_vec[ri] = min_dists[idx]
        cluster_sizes[c] += 1

        key = (c, cdr3)
        cnt = get(cluster_cdr3_counts, key, 0) + 1
        cluster_cdr3_counts[key] = cnt
        if cnt > cluster_max_cdr3[c]
            cluster_max_cdr3[c] = cnt
        end
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
    collapse_lineages(df::DataFrame, strategy=Hardest(); tie_breaker, tie_atol, abundance) -> DataFrame

Collapse lineages to representative sequences.

# Arguments
- `df::DataFrame`: Data from `process_lineages` with `lineage_id` column
- `strategy::CollapseStrategy=Hardest()`: Collapse strategy
  - `Hardest()`: One representative per lineage (highest clone frequency, then `tie_breaker`)
  - `IgDiscover()`: IgDiscover `clonotypes` representative (lineage-wide most common `vdj_nt`)
  - `Soft(cutoff)`: Keep all clones with frequency ≥ cutoff

# Keyword Arguments
- `tie_breaker::AbstractTieBreaker`: Strategy for breaking ties. Default depends on `strategy`.
- `tie_atol::Real=0.0`: Tolerance for frequency comparison (`Hardest` only)
- `abundance::Symbol=:count`: How to interpret abundance.
  - `:count` — use the `count` column when present (missing treated as 1); otherwise each row is 1.
    Values are treated as **additive** across lineage members (read/UMI/observation counts).
  - `:rows` — ignore `count`; each row contributes 1.

Clone frequency used by `Hardest`/`Soft` is the fraction of **rows** in a clone, not the
sum of `count`.

# Returns
For `Hardest()` and `IgDiscover()`:
- Selected rows with added `count::Int` (sum of abundances) and `nVDJ_nt::Int` (unique VDJ sequences)

For `Soft(cutoff)`:
- Rows meeting threshold with added `clone_frequency::Float64` and `sequence_count::Int`
"""
function collapse_lineages(df::DataFrame, strategy::CollapseStrategy=Hardest();
                           tie_breaker::AbstractTieBreaker=default_tie_breaker(strategy),
                           tie_atol::Real=0.0,
                           abundance::Symbol=:count)
    abundance = validate_abundance(abundance)
    nrow(df) == 0 && return empty_collapse_result(df, strategy)
    validate_tie_breaker(tie_breaker, df)
    collapse_lineages_impl(df, strategy, tie_breaker, Float64(tie_atol), abundance)
end

function empty_collapse_result(df::DataFrame, ::Hardest)
    result = copy(df)
    result[!, :count] = Int[]
    result[!, :nVDJ_nt] = Int[]
    result
end

function empty_collapse_result(df::DataFrame, ::IgDiscover)
    empty_collapse_result(df, Hardest())
end

function empty_collapse_result(df::DataFrame, ::Soft)
    result = copy(df)
    result[!, :clone_frequency] = Float64[]
    result[!, :sequence_count] = Int[]
    result
end

function collapse_lineages_impl(df::DataFrame, strategy::Hardest, tie_breaker::AbstractTieBreaker,
                                 atol::Float64, abundance::Symbol)
    clone_lookup = compute_clone_frequency_lookup(df)
    vdj_lookup = compute_vdj_count_lookup(tie_breaker, df, abundance)
    agg = compute_aggregates(df, abundance)

    groups = groupby(df, :lineage_id; sort=false)
    selected = Vector{Int}(undef, length(groups))
    for (gi, group) in enumerate(groups)
        rows = parentindices(group)[1]
        selected[gi] = select_hardest_index(df, rows, tie_breaker, atol, clone_lookup, vdj_lookup, abundance)
    end

    build_hardest_result(df, selected, agg)
end

function collapse_lineages_impl(df::DataFrame, ::IgDiscover, tie_breaker::AbstractTieBreaker,
                                 ::Float64, abundance::Symbol)
    vdj_lookup = NoVdjCountLookup()
    clone_lookup = CloneFrequencyLookup(Dict{CloneKey,Float64}(), Dict{CloneKey,Int}())
    agg = compute_aggregates(df, abundance)

    groups = groupby(df, :lineage_id; sort=false)
    selected = Vector{Int}(undef, length(groups))
    for (gi, group) in enumerate(groups)
        rows = parentindices(group)[1]
        selected[gi] = select_representative(tie_breaker, df, rows, collect(Int, rows),
                                           clone_lookup, vdj_lookup, abundance)
    end

    build_hardest_result(df, selected, agg)
end

function collapse_lineages_impl(df::DataFrame, strategy::Soft, ::AbstractTieBreaker, ::Float64, ::Symbol)
    clone_lookup = compute_clone_frequency_lookup(df)
    cutoff = strategy.cutoff

    indices = Int[]
    sizehint!(indices, nrow(df) ÷ 2)
    @inbounds for i in 1:nrow(df)
        key = (df.d_region[i], df.lineage_id[i], df.j_call_first[i], df.v_call_first[i], df.cdr3[i])::CloneKey
        clone_lookup.freq[key] >= cutoff && push!(indices, i)
    end

    build_soft_result(df, indices, clone_lookup)
end

# Validation via dispatch
validate_tie_breaker(::MostCommonVdjNtTieBreaker, df::DataFrame) =
    hasproperty(df, :vdj_nt) || throw(ArgumentError("vdj_nt column required"))

validate_tie_breaker(::LineageMostCommonVdjNtTieBreaker, df::DataFrame) =
    hasproperty(df, :vdj_nt) || throw(ArgumentError("vdj_nt column required"))

validate_tie_breaker(::LowestVShmTieBreaker, df::DataFrame) =
    hasproperty(df, :V_SHM) || throw(ArgumentError("V_SHM column required for ByLowestVShm()"))

function validate_tie_breaker(tb::TieBreaker, df::DataFrame)
    any(col == :vdj_count for (col, _) in tb.criteria) && !hasproperty(df, :vdj_nt) &&
        throw(ArgumentError("vdj_nt column required for vdj_count"))
end

validate_tie_breaker(::AbstractTieBreaker, ::DataFrame) = nothing

# VDJ count lookup via dispatch
compute_vdj_count_lookup(::MostCommonVdjNtTieBreaker, ::DataFrame, ::Symbol) = NoVdjCountLookup()
compute_vdj_count_lookup(::LineageMostCommonVdjNtTieBreaker, ::DataFrame, ::Symbol) = NoVdjCountLookup()
compute_vdj_count_lookup(::LowestVShmTieBreaker, ::DataFrame, ::Symbol) = NoVdjCountLookup()

function compute_vdj_count_lookup(tb::TieBreaker, df::DataFrame, abundance::Symbol)
    any(col == :vdj_count for (col, _) in tb.criteria) ? compute_vdj_count_lookup_impl(df, abundance) : NoVdjCountLookup()
end

function compute_clone_frequency_lookup(df::DataFrame)
    clone_keys = [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3]
    clone_counts = combine(groupby(df, clone_keys; sort=false), nrow => :_cnt)
    lineage_totals = combine(groupby(clone_counts, :lineage_id; sort=false), :_cnt => sum => :_total)

    isempty(lineage_totals) && return CloneFrequencyLookup(Dict{CloneKey,Float64}(), Dict{CloneKey,Int}())

    n_lineages = maximum(lineage_totals.lineage_id)
    lin_totals = zeros(Int, n_lineages)
    @inbounds for row in eachrow(lineage_totals)
        lin_totals[row.lineage_id] = row._total
    end

    freq = Dict{CloneKey,Float64}()
    cnt_dict = Dict{CloneKey,Int}()

    for row in eachrow(clone_counts)
        key = (row.d_region, row.lineage_id, row.j_call_first, row.v_call_first, row.cdr3)::CloneKey
        freq[key] = row._cnt / lin_totals[row.lineage_id]
        cnt_dict[key] = row._cnt
    end

    CloneFrequencyLookup(freq, cnt_dict)
end

function compute_vdj_count_lookup_impl(df::DataFrame, abundance::Symbol)
    grouped = combine(groupby(df, [:lineage_id, :d_region, :j_call_first, :v_call_first, :cdr3, :vdj_nt]; sort=false)) do g
        s = 0
        for i in 1:nrow(g)
            s += row_abundance(g, i, abundance)
        end
        (; _vdj_cnt = s)
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

function select_hardest_index(df::DataFrame, rows::AbstractVector{Int},
                              tie_breaker::AbstractTieBreaker, atol::Float64,
                              clone_lookup::CloneFrequencyLookup, vdj_lookup, abundance::Symbol)
    if !uses_frequency_candidates(tie_breaker)
        return select_representative(tie_breaker, df, rows, collect(Int, rows), clone_lookup, vdj_lookup, abundance)
    end

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

    select_representative(tie_breaker, df, rows, candidates, clone_lookup, vdj_lookup, abundance)
end

function select_by_vdj_abundance(df::DataFrame, rows::AbstractVector{Int}, abundance::Symbol;
                                small_lineage_first_row::Bool, prefer_lex_vdj::Bool)
    n = length(rows)
    small_lineage_first_row && n <= 2 && return rows[1]

    vdj_col = df.vdj_nt
    vdj_counts = Dict{String,Int}()
    vdj_first_row = Dict{String,Int}()

    for ri in rows
        @inbounds vdj = vdj_col[ri]
        ismissing(vdj) && continue
        cnt = row_abundance(df, ri, abundance)

        if haskey(vdj_counts, vdj)
            vdj_counts[vdj] += cnt
        else
            vdj_counts[vdj] = cnt
            vdj_first_row[vdj] = ri
        end
    end

    isempty(vdj_counts) && return rows[1]

    max_count = -1
    best_row = rows[1]
    best_vdj = ""

    for (vdj, cnt) in vdj_counts
        first_row = vdj_first_row[vdj]
        better = if cnt != max_count
            cnt > max_count
        elseif prefer_lex_vdj
            vdj < best_vdj
        else
            first_row < best_row
        end
        if better
            max_count = cnt
            best_row = first_row
            best_vdj = vdj
        end
    end

    best_row
end

# Representative selection via dispatch on tie_breaker type
function select_representative(::MostCommonVdjNtTieBreaker, df::DataFrame,
                              ::AbstractVector{Int}, candidates::Vector{Int},
                              ::CloneFrequencyLookup, ::NoVdjCountLookup, abundance::Symbol)
    select_by_vdj_abundance(df, candidates, abundance;
                           small_lineage_first_row=false, prefer_lex_vdj=true)
end

function select_representative(tb::LineageMostCommonVdjNtTieBreaker, df::DataFrame,
                              rows::AbstractVector{Int}, ::Vector{Int},
                              ::CloneFrequencyLookup, ::NoVdjCountLookup, abundance::Symbol)
    select_by_vdj_abundance(df, rows, abundance;
                             small_lineage_first_row=tb.small_lineage_first_row,
                             prefer_lex_vdj=false)
end

function select_representative(::LowestVShmTieBreaker, df::DataFrame,
                              rows::AbstractVector{Int}, ::Vector{Int},
                              ::CloneFrequencyLookup, ::NoVdjCountLookup, abundance::Symbol)
    shm = df.V_SHM
    has_vdj = hasproperty(df, :vdj_nt)
    best = rows[1]
    @inbounds for k in 2:length(rows)
        idx = rows[k]
        if is_better_vshm(df, idx, best, shm, has_vdj, abundance)
            best = idx
        end
    end
    best
end

function is_better_vshm(df, idx1::Int, idx2::Int, shm, has_vdj::Bool, abundance::Symbol)
    s1 = shm[idx1]
    s2 = shm[idx2]
    m1 = ismissing(s1)
    m2 = ismissing(s2)
    if m1 && !m2
        return false
    elseif m2 && !m1
        return true
    elseif !m1 && !m2 && s1 != s2
        return s1 < s2
    end

    c1 = row_abundance(df, idx1, abundance)
    c2 = row_abundance(df, idx2, abundance)
    c1 != c2 && return c1 > c2

    if has_vdj
        v1 = df.vdj_nt[idx1]
        v2 = df.vdj_nt[idx2]
        if !isequal(v1, v2)
            ismissing(v1) && return false
            ismissing(v2) && return true
            return v1 < v2
        end
    end

    c3a = df.cdr3[idx1]
    c3b = df.cdr3[idx2]
    isequal(c3a, c3b) && return false
    return c3a < c3b
end

function select_representative(tb::TieBreaker, df::DataFrame,
                              ::AbstractVector{Int}, candidates::Vector{Int},
                              clone_lookup::CloneFrequencyLookup, vdj_lookup, ::Symbol)
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

# Aggregation for Hardest collapse (uses Vectors since lineage_id is 1:N)
struct LineageAggregates
    count_sums::Vector{Int}
    vdj_counts::Vector{Int}
end

function compute_aggregates(df::AbstractDataFrame, abundance::Symbol=:count)
    has_vdj = hasproperty(df, :vdj_nt)
    nrow(df) == 0 && return LineageAggregates(Int[], Int[])

    n_lineages = maximum(df.lineage_id)
    count_sums = zeros(Int, n_lineages)
    vdj_sets = [Set{String}() for _ in 1:n_lineages]

    lineage_col = df.lineage_id

    @inbounds for i in 1:nrow(df)
        lid = lineage_col[i]
        count_sums[lid] += row_abundance(df, i, abundance)

        if has_vdj
            vdj = df.vdj_nt[i]
            !ismissing(vdj) && push!(vdj_sets[lid], vdj)
        end
    end

    vdj_counts = has_vdj ? [length(s) for s in vdj_sets] : Int[]
    LineageAggregates(count_sums, vdj_counts)
end

function build_hardest_result(df::DataFrame, selected::Vector{Int}, agg::LineageAggregates)
    result = df[selected, :]
    n = length(selected)
    n == 0 && return empty_collapse_result(df, Hardest())

    counts = Vector{Int}(undef, n)
    has_vdj = !isempty(agg.vdj_counts)
    nvdj = has_vdj ? Vector{Int}(undef, n) : Int[]

    @inbounds for i in 1:n
        lid = result.lineage_id[i]
        counts[i] = agg.count_sums[lid]
        has_vdj && (nvdj[i] = agg.vdj_counts[lid])
    end

    result[!, :count] = counts
    result[!, :nVDJ_nt] = has_vdj ? nvdj : fill(missing, n)
    result
end

function build_soft_result(df::DataFrame, selected::Vector{Int}, clone_lookup::CloneFrequencyLookup)
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
    nrow(df) == 0 && return DataFrame(
        lineage_id=Int[], max_clone_frequency=Float64[],
        hardest_tie_count=Int[], hardest_tied=Bool[], hardest_tied_cdr3s=Vector{Vector{String}}()
    )
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
