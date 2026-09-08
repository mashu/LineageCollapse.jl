"""
    extract_interval(seq, start_idx, stop_idx)

Return `seq[start_idx:stop_idx]`, or `missing` when coordinates are absent or out of bounds.
Empty intervals (`start_idx > stop_idx`) are invalid for VDJ spans.
"""
function extract_interval(seq, start_idx, stop_idx)
    (ismissing(seq) || ismissing(start_idx) || ismissing(stop_idx)) && return missing
    n = length(seq)
    start_i = Int(start_idx)
    stop_i = Int(stop_idx)
    (1 <= start_i <= stop_i <= n) || return missing
    return seq[start_i:stop_i]
end

"""
    extract_d_region(seq, v_end, j_start)

D-region slice `seq[v_end+1:j_start]`. `v_end` and `j_start` must be in-bounds
1-based AIRR coordinates. An empty D (V immediately followed by J) is allowed.
"""
function extract_d_region(seq, v_end, j_start)
    (ismissing(seq) || ismissing(v_end) || ismissing(j_start)) && return missing
    n = length(seq)
    v = Int(v_end)
    j = Int(j_start)
    (1 <= v <= n && 1 <= j <= n) || return missing
    start_i = v + 1
    start_i > j && return ""
    return seq[start_i:j]
end

"""
    first_allele(call)

Return the first allele token from an AIRR gene call.

Comma-separated multiple assignments (`IGHV1-2*01,IGHV1-2*02`) keep the first
token. Slash-joined dual locus names (`IGHV3-23*01/IGHV3-23D*01`) are kept
intact because they are a single assignment, not two alternatives.
"""
function first_allele(call)
    ismissing(call) && return missing
    return String(first(split(string(call), ",")))
end

"""
    deduplicate_data(df::DataFrame, use_barcode::Bool=false)::DataFrame

Deduplicate the input DataFrame based on sequence or sequence+barcode.

# Arguments
- `df::DataFrame`: Input DataFrame.
- `use_barcode::Bool=false`: Whether to use barcode for deduplication.

# Returns
- `DataFrame`: Deduplicated DataFrame.
"""
function deduplicate_data(df::DataFrame, use_barcode::Bool=false)::DataFrame
    if use_barcode && :barcode in propertynames(df)
        df = unique(df, [:sequence, :barcode])
        @info "Deduplicated using sequence and barcode: $(nrow(df))"
    else
        df = unique(df, :sequence)
        @info "Deduplicated using sequence only: $(nrow(df))"
    end
    return df
end

"""
    preprocess_data(df::DataFrame; min_d_region_length=nothing, deduplicate=false, use_barcode=false, rename_count_to_consensus_reads=false)::DataFrame

Preprocess the input DataFrame by performing data cleaning and transformation.

# Arguments
- `df::DataFrame`: Input DataFrame.
- `min_d_region_length::Union{Int,Nothing}=nothing`: Inclusive minimum D-region length.
  Rows with `length(d_region) < min_d_region_length` are dropped. If `nothing`, no D-length filter.
- `deduplicate::Bool=false`: Whether to deduplicate the DataFrame.
- `use_barcode::Bool=false`: Whether to use barcode for deduplication (only applicable if deduplicate is true).
- `rename_count_to_consensus_reads::Bool=false`: If true, renames `:count` to `:consensus_reads` before downstream steps. Use this when `:count` holds PCR- or consensus-level abundances that must not be summed like per-sequence read counts during lineage collapse.

If the input has a `:count` column with any value greater than 1 and `rename_count_to_consensus_reads` is false, a warning is emitted: lineage collapse aggregates `:count` across members, which can undo PCR bias correction if `:count` was already consensus-level abundance.

Rows with an explicit `stop_codon == true` are dropped. Missing `stop_codon` is treated as
*not* having a stop codon (the row is kept). Rows with missing `cdr3`, `v_call`, or `j_call`
are dropped. Invalid V/J coordinates yield `missing` D-region or `vdj_nt` and those rows
are dropped.

An empty input yields an empty DataFrame with the derived columns present.
"""
function preprocess_data(df::DataFrame;
    min_d_region_length::Union{Int,Nothing}=nothing,
    deduplicate::Bool=false,
    use_barcode::Bool=false,
    rename_count_to_consensus_reads::Bool=false)::DataFrame
    if hasproperty(df, :count)
        if rename_count_to_consensus_reads
            rename!(df, :count => :consensus_reads)
        else
            vals = collect(skipmissing(df.count))
            if !isempty(vals) && maximum(vals) > 1
                @warn "Input has `:count` with values > 1. During lineage collapse, `count` is summed across lineage members (abundance=:count). If this column holds consensus-level or PCR-corrected abundances (not additive sequence counts), rename it before processing: pass `rename_count_to_consensus_reads=true` to `preprocess_data`, or call `rename!(df, :count => :consensus_reads)` so aggregation does not undo bias correction. Use `collapse_lineages(..., abundance=:rows)` to ignore `count`."
            end
        end
    end

    @info "Processing $(nrow(df)) rows"

    nrow(df) == 0 && return empty_preprocessed(df)

    # Remove rows with missing CDR3
    df = dropmissing(df, :cdr3)
    @info "Dropped missing CDR3 rows: $(nrow(df))"

    # Drop explicit stop codons; keep missing productivity annotations
    if hasproperty(df, :stop_codon)
        df = filter(row -> coalesce(row.stop_codon, false) == false, df)
        @info "Dropped stop codons: $(nrow(df))"
    end

    nrow(df) == 0 && return empty_preprocessed(df)

    # Calculate D region with explicit coordinate validation
    transform!(df,
        [:sequence, :v_sequence_end, :j_sequence_start] =>
        ByRow(extract_d_region) => :d_region
    )
    df = dropmissing(df, :d_region)

    # Derive VDJ_nt if missing and coordinates are available (MIAIRR input)
    if !(:vdj_nt in propertynames(df)) &&
       (:v_sequence_start in propertynames(df)) &&
       (:j_sequence_end in propertynames(df)) &&
       (:sequence in propertynames(df))
        transform!(df,
            [:sequence, :v_sequence_start, :j_sequence_end] =>
            ByRow((seq, v_start, j_end) -> extract_interval(seq, v_start, j_end)) => :vdj_nt
        )
    end

    # Filter based on D region length if specified (inclusive minimum)
    if !isnothing(min_d_region_length)
        df = filter(row -> length(row.d_region) >= min_d_region_length, df)
        @info "Dropped D region sequences shorter than $min_d_region_length: $(nrow(df))"
    end

    # Extract first allele from v_call and j_call (normalize to String)
    transform!(df, :v_call => ByRow(first_allele) => :v_call_first)
    transform!(df, :j_call => ByRow(first_allele) => :j_call_first)
    df = dropmissing(df, [:v_call_first, :j_call_first])

    nrow(df) == 0 && return df

    # Normalize string columns to String type for consistent downstream processing
    df[!, :d_region] = String.(df.d_region)
    df[!, :cdr3] = String.(df.cdr3)

    # Calculate CDR3 length
    transform!(df, :cdr3 => ByRow(length) => :cdr3_length)

    # Deduplicate if requested
    if deduplicate
        df = deduplicate_data(df, use_barcode)
    end

    return df
end

function empty_preprocessed(df::DataFrame)
    result = copy(df)
    :d_region in propertynames(result) || (result[!, :d_region] = String[])
    :v_call_first in propertynames(result) || (result[!, :v_call_first] = String[])
    :j_call_first in propertynames(result) || (result[!, :j_call_first] = String[])
    :cdr3_length in propertynames(result) || (result[!, :cdr3_length] = Int[])
    return result
end
