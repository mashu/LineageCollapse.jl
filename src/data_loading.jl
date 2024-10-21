"""
    load_data(filepath::String; 
              delimiter::Char='\t', 
              required_columns=[:sequence_id, :sequence, :v_sequence_end, :j_sequence_start, :cdr3, :v_call, :j_call, :stop_codon])::DataFrame

Load data from a file (compressed or uncompressed) and return a DataFrame.

# Arguments
- `filepath::String`: Path to the data file.
- `delimiter::Char='\t'`: Delimiter used in the data file (default: tab).
- `required_columns::Vector{Symbol}`: Required columns to select from the data file.

# Returns
- `DataFrame`: DataFrame containing the loaded data.

# Throws
- `ArgumentError`: If any of the required columns are missing in the data file.
"""
function load_data(filepath::String; 
                   delimiter::Char='\t', 
                   required_columns=[:sequence_id, :sequence, :v_sequence_end, :j_sequence_start, :cdr3, :v_call, :j_call, :stop_codon])::DataFrame
    try
        # Read all columns
        df = CSV.File(filepath, delim=delimiter) |> DataFrame
        available_columns = names(df)

        # Check if required columns are present
        missing_columns = setdiff(required_columns, Symbol.(available_columns))
        if !isempty(missing_columns)
            throw(ArgumentError("Missing required columns: $(join(missing_columns, ", "))"))
        end

        return df
    catch e
        throw(ErrorException("Error loading data: $(sprint(showerror, e))"))
    end
end