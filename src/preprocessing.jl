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
    preprocess_data(df::DataFrame; min_d_region_length::Union{Int,Nothing}=nothing, deduplicate::Bool=false, use_barcode::Bool=false, rename_count_to_consensus_reads::Bool=false)::DataFrame

Preprocess the input DataFrame by performing data cleaning and transformation.

# Arguments
- `df::DataFrame`: Input DataFrame.
- `min_d_region_length::Union{Int,Nothing}=nothing`: Minimum length of the D region to keep. If nothing, no filtering is applied.
- `deduplicate::Bool=false`: Whether to deduplicate the DataFrame.
- `use_barcode::Bool=false`: Whether to use barcode for deduplication (only applicable if deduplicate is true).
- `rename_count_to_consensus_reads::Bool=false`: If true, renames `:count` to `:consensus_reads` before downstream steps. Use this when `:count` holds PCR- or consensus-level abundances that must not be summed like per-sequence read counts during lineage collapse.

If the input has a `:count` column with any value greater than 1 and `rename_count_to_consensus_reads` is false, a warning is emitted: lineage collapse aggregates `:count` across members, which can undo PCR bias correction if `:count` was already consensus-level abundance.

# Returns
- `DataFrame`: Preprocessed DataFrame.
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
                @warn "Input has `:count` with values > 1. During lineage collapse, `count` is summed across lineage members. If this column holds consensus-level or PCR-corrected abundances (not additive sequence counts), rename it before processing: pass `rename_count_to_consensus_reads=true` to `preprocess_data`, or call `rename!(df, :count => :consensus_reads)` so aggregation does not undo bias correction."
            end
        end
    end

    @info "Processing $(nrow(df)) rows"
    
    # Remove rows with missing CDR3
    df = dropmissing(df, :cdr3)
    @info "Dropped missing CDR3 rows: $(nrow(df))"
    
    # Remove rows with stop codons
    df = filter(row -> row.stop_codon == false, df)
    @info "Dropped stop codons: $(nrow(df))"
    
    # Calculate D region
    transform!(df,
        [:sequence, :v_sequence_end, :j_sequence_start] =>
        ByRow((seq, v_end, j_start) -> seq[v_end+1:j_start]) => :d_region
    )

    # Derive VDJ_nt if missing and coordinates are available (MIAIRR input)
    if !(:vdj_nt in propertynames(df)) &&
       (:v_sequence_start in propertynames(df)) &&
       (:j_sequence_end in propertynames(df)) &&
       (:sequence in propertynames(df))
        transform!(df,
            [:sequence, :v_sequence_start, :j_sequence_end] =>
            ByRow((seq, v_start, j_end) -> begin
                if ismissing(seq) || ismissing(v_start) || ismissing(j_end)
                    return missing
                end
                if !(1 <= v_start <= j_end <= length(seq))
                    return missing
                end
                return seq[v_start:j_end]
            end) => :vdj_nt
        )
    end

    # Filter based on D region length if specified
    if !isnothing(min_d_region_length)
        df = filter(row -> length(row.d_region) > min_d_region_length, df)
        @info "Dropped short (≤$min_d_region_length) D region sequences: $(nrow(df))"
    end
    
    # Extract first allele from v_call and j_call (normalize to String)
    transform!(df, :v_call => ByRow(x -> String(first(split(x, ",")))) => :v_call_first)
    transform!(df, :j_call => ByRow(x -> String(first(split(x, ",")))) => :j_call_first)
    
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
