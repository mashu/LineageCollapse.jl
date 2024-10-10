"""
    preprocess_data(df::DataFrame; min_d_region_length::Int=9)::DataFrame

Preprocess the input DataFrame by performing data cleaning and transformation.

# Arguments
- `df::DataFrame`: Input DataFrame.
- `min_d_region_length::Int=9`: Minimum length of the D sequence to keep.

# Returns
- `DataFrame`: Preprocessed DataFrame.
"""
function preprocess_data(df::DataFrame; min_d_region_length::Int=9)::DataFrame
    @info "Processing $(nrow(df)) rows"
    dropmissing!(df, [:cdr3])
    @info "Dropped missing CDR3 rows: $(nrow(df))"
    filter!(x->x.stop_codon == false, df)
    @info "Dropped stop codons: $(nrow(df))"
    unique!(df, :sequence)
    @info "Dropped duplicated sequences: $(nrow(df))"

    transform!(df,
        [:sequence, :v_sequence_end, :j_sequence_start] =>
        ByRow((seq, v_end, j_start) -> seq[v_end:j_start]) => :d_region,
    )

    filter!(row -> length(row.d_region) > min_d_region_length, df)
    @info "Dropped short (<$min_d_region_length) D region sequences: $(nrow(df))"
    
    df = transform(df, :v_call => ByRow(x -> first(split(x, ","))) => :v_call_first)
    df = transform(df, :j_call => ByRow(x -> first(split(x, ","))) => :j_call_first)
    transform!(df, :cdr3 => ByRow(length) => :cdr3_length)

    return df
end
