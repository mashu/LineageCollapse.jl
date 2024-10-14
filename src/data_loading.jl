"""
    load_data(filepath::String; 
              delimiter::Char='\t', 
              required_columns=[:sequence_id, :sequence, :v_sequence_end, :j_sequence_start, :cdr3, :v_call, :j_call, :stop_codon])::DataFrame

Load data from a file and return a DataFrame.

# Arguments
- `filepath::String`: Path to the data file.
- `delimiter::Char='\t'`: Delimiter used in the data file (default: tab).
- `required_columns::Vector{Symbol}`: Required columns in the data file.

# Returns
- `DataFrame`: DataFrame containing the loaded data.

# Throws
- `ArgumentError`: If any of the required columns are missing in the data file.
"""
function load_data(filepath::String; 
                   delimiter::Char='\t', 
                   required_columns=[:sequence_id, :sequence, :v_sequence_end, :j_sequence_start, :cdr3, :v_call, :j_call, :stop_codon])::DataFrame
    # First, read the header to check for required columns
    header = String.(split(readline(filepath), delimiter))
    missing_columns = setdiff(String.(required_columns), header)
    
    if !isempty(missing_columns)
        throw(ArgumentError("Missing required columns: $(join(missing_columns, ", "))"))
    end

    # If all required columns are present, load the data
    return CSV.File(filepath, delim=delimiter, select=required_columns) |> DataFrame
end