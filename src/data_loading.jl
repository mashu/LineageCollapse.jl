"""
    load_data(filepath::String; delimiter::Char='\t', required_columns=[:sequence, :v_sequence_end, :j_sequence_start, :cdr3, :v_call, :j_call, :stop_codon])::DataFrame

Load data from a file and return a DataFrame.

# Arguments
- `filepath::String`: Path to the data file.
- `delimiter::Char='\t'`: Delimiter used in the data file (default: tab).
- `required_columns::Vector{Symbol}=[:sequence, :v_sequence_end, :j_sequence_start, :cdr3, :v_call, :j_call, :stop_codon]`: Required columns in the data file.

# Returns
- `DataFrame`: DataFrame containing the loaded data.
"""
function load_data(filepath::String; 
                   delimiter::Char='\t', 
                   required_columns=[:sequence, :v_sequence_end, :j_sequence_start, :cdr3, :v_call, :j_call, :stop_codon])::DataFrame
    return CSV.File(filepath, delim=delimiter, select=required_columns) |> DataFrame
end
