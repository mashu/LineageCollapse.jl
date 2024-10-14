"""
    preprocess_data(df::DataFrame; min_d_region_length::Int=0)::DataFrame

Preprocess the input DataFrame by performing data cleaning and transformation.

# Arguments
- `df::DataFrame`: Input DataFrame.
- `min_d_region_length::Int=0`: Minimum length of the D region to keep.

# Returns
- `DataFrame`: Preprocessed DataFrame.
"""
function preprocess_data(df::DataFrame; min_d_region_length::Int=0)::DataFrame
    @info "Processing $(nrow(df)) rows"
    
    # Remove rows with missing CDR3
    df = dropmissing(df, :cdr3)
    @info "Dropped missing CDR3 rows: $(nrow(df))"
    
    # Remove rows with stop codons
    df = filter(row -> row.stop_codon == false, df)
    @info "Dropped stop codons: $(nrow(df))"
    
    # Remove duplicate sequences
    df = unique(df, :sequence)
    @info "Dropped duplicated sequences: $(nrow(df))"

    # Calculate D region
    transform!(df,
        [:sequence, :v_sequence_end, :j_sequence_start] =>
        ByRow((seq, v_end, j_start) -> seq[v_end+1:j_start]) => :d_region
    )

    # Filter based on D region length
    df = filter(row -> length(row.d_region) > min_d_region_length, df)
    @info "Dropped short (≤$min_d_region_length) D region sequences: $(nrow(df))"
    
    # Extract first allele from v_call and j_call
    transform!(df, :v_call => ByRow(x -> first(split(x, ","))) => :v_call_first)
    transform!(df, :j_call => ByRow(x -> first(split(x, ","))) => :j_call_first)
    
    # Calculate CDR3 length
    transform!(df, :cdr3 => ByRow(length) => :cdr3_length)

    return df
end