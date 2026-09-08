module LineageCollapseCLI

using ArgParse
using CSV
using DataFrames
using LineageCollapse
using Logging

export main, julia_main, run_pipeline

abstract type AbstractFilterOperator end
struct LessThan <: AbstractFilterOperator end
struct LessEqual <: AbstractFilterOperator end
struct GreaterThan <: AbstractFilterOperator end
struct GreaterEqual <: AbstractFilterOperator end
struct EqualTo <: AbstractFilterOperator end
struct NotEqualTo <: AbstractFilterOperator end

struct ColumnFilter{O<:AbstractFilterOperator,T}
    column::Symbol
    value::T
end

struct CliConfig{T<:Real}
    input::String
    output::String
    delimiter::Char
    default_filters::Bool
    filters::Vector{String}
    rename_count_to::Union{Nothing,Symbol}
    min_d_region_length::Union{Nothing,Int}
    deduplicate::Bool
    use_barcode::Bool
    threshold::T
    linkage::Symbol
    collapse_strategy::String
    soft_cutoff::Float64
    tie_breaker::String
    tie_atol::Float64
    abundance::Symbol
    processed_output::Union{Nothing,String}
    quiet::Bool
end

const DEFAULT_FILTER_SPECS = ["V_errors<=1", "J_errors<=0"]
const FILTER_OPERATORS = [
    "<=" => LessEqual,
    ">=" => GreaterEqual,
    "==" => EqualTo,
    "!=" => NotEqualTo,
    "<" => LessThan,
    ">" => GreaterThan,
]
const LINKAGES = Set([:single, :complete, :average])
const COLLAPSE_STRATEGIES = Set(["hardest", "soft", "igdiscover"])
const TIE_BREAKERS = Set([
    "most-common-vdj-nt",
    "lineage-most-common-vdj-nt",
    "lowest-v-shm",
    "vdj-count",
    "cdr3-count",
    "sequence-count",
    "lexicographic",
    "first",
    "most-naive",
])

function build_parser()
    settings = ArgParseSettings(
        prog = "lineage-collapse",
        description = "Read one TSV file, filter, process CDR3 lineages, and write collapsed lineage output.",
    )

    @add_arg_table! settings begin
        "input"
            help = "Input TSV/CSV file path. Compressed files supported by CSV.jl."
            required = true
            arg_type = String
        "--output", "-o"
            help = "Output path for collapsed lineages. Use '-' for stdout."
            default = "-"
        "--processed-output"
            help = "Optional path for the processed pre-collapse lineage table."
            default = nothing
        "--delimiter"
            help = "Input/output delimiter. Use '\\t' for tab."
            default = "\\t"
        "--no-default-filters"
            help = "Disable default filters V_errors<=1 and J_errors<=0."
            action = :store_true
        "--filter"
            help = "Additional row filter such as 'productive==true', 'V_errors<=1', or 'cdr3_length>8'. May be repeated."
            action = :append_arg
            default = String[]
        "--keep-count"
            help = "Keep an existing count column as additive sequence counts. By default count is renamed."
            action = :store_true
        "--rename-count-to"
            help = "Column name used when renaming an existing count column."
            default = "consensus_count"
        "--min-d-region-length"
            help = "Drop rows whose inferred D region length is strictly less than this value (inclusive minimum)."
            arg_type = Int
            default = nothing
        "--deduplicate"
            help = "Deduplicate by sequence before lineage processing."
            action = :store_true
        "--use-barcode"
            help = "When deduplicating, use sequence and barcode if barcode exists."
            action = :store_true
        "--cdr3-differences"
            help = "Maximum absolute CDR3 differences for clustering."
            arg_type = Int
            default = 1
        "--cdr3-distance-fraction"
            help = "Maximum normalized CDR3 distance, from 0.0 to 1.0. Overrides --cdr3-differences."
            arg_type = Float64
            default = nothing
        "--cdr3-similarity-percent"
            help = "Minimum CDR3 similarity percent, from 0 to 100. Overrides --cdr3-differences."
            arg_type = Float64
            default = nothing
        "--linkage"
            help = "Hierarchical clustering linkage: single, complete, or average."
            default = "single"
        "--collapse"
            help = "Collapse output strategy: hardest, soft, or igdiscover."
            default = "hardest"
        "--soft-cutoff"
            help = "For --collapse soft, keep clones with frequency at least this value."
            arg_type = Float64
            default = 0.3
        "--tie-breaker"
            help = """
Hardest representative policy (default: most-common-vdj-nt).
  most-common-vdj-nt — among maximum clone-frequency clones, pick the vdj_nt with largest summed abundance.
  lineage-most-common-vdj-nt — IgDiscover-style: sum count for each vdj_nt across the whole lineage (not only top-frequency clones).
  lowest-v-shm — least-mutated V (minimum V_SHM) among all lineage members.
  vdj-count, cdr3-count, sequence-count, lexicographic, first, most-naive — see package docs."""
            default = "most-common-vdj-nt"
        "--tie-atol"
            help = "Absolute tolerance when comparing clone frequencies for Hardest ties."
            arg_type = Float64
            default = 0.0
        "--abundance"
            help = "How to interpret abundance: count (use the count column when present; missing=1) or rows (each row is 1)."
            default = "count"
        "--quiet"
            help = "Suppress progress messages."
            action = :store_true
    end

    return settings
end

function delimiter_char(delimiter::AbstractString)
    delimiter == "\\t" && return '\t'
    length(delimiter) == 1 || throw(ArgumentError("--delimiter must be one character or '\\t'"))
    return only(delimiter)
end

function selected_threshold(args)
    similarity = args["cdr3-similarity-percent"]
    distance = args["cdr3-distance-fraction"]
    !isnothing(similarity) && !isnothing(distance) &&
        throw(ArgumentError("Use either --cdr3-similarity-percent or --cdr3-distance-fraction, not both"))

    if !isnothing(similarity)
        0.0 <= similarity <= 100.0 ||
            throw(ArgumentError("--cdr3-similarity-percent must be between 0 and 100"))
        return 1.0 - similarity / 100.0
    end

    if !isnothing(distance)
        0.0 <= distance <= 1.0 ||
            throw(ArgumentError("--cdr3-distance-fraction must be between 0.0 and 1.0"))
        return distance
    end

    differences = args["cdr3-differences"]
    differences >= 0 || throw(ArgumentError("--cdr3-differences must be non-negative"))
    return differences
end

function config_from_args(args)
    linkage = Symbol(String(args["linkage"]))
    linkage in LINKAGES || throw(ArgumentError("--linkage must be one of single, complete, or average"))

    collapse_strategy = String(args["collapse"])
    collapse_strategy in COLLAPSE_STRATEGIES || throw(ArgumentError("--collapse must be hardest, soft, or igdiscover"))

    tie_breaker = String(args["tie-breaker"])
    tie_breaker in TIE_BREAKERS || throw(ArgumentError("Unknown --tie-breaker: $tie_breaker"))

    soft_cutoff = Float64(args["soft-cutoff"])
    0.0 <= soft_cutoff <= 1.0 || throw(ArgumentError("--soft-cutoff must be between 0.0 and 1.0"))

    rename_count_to = args["keep-count"] ? nothing : Symbol(String(args["rename-count-to"]))

    input = args["input"]
    (input isa AbstractString && !isempty(strip(input))) ||
        throw(ArgumentError("input path must be a non-empty string"))

    filters_raw = get(args, "filter", String[])
    filters = filters_raw isa Vector{String} ? filters_raw : String[String(x) for x in filters_raw]

    proc = get(args, "processed-output", nothing)
    processed_output = if isnothing(proc) || (proc isa AbstractString && isempty(strip(proc)))
        nothing
    else
        String(proc)
    end

    abundance = Symbol(String(args["abundance"]))
    (abundance === :count || abundance === :rows) ||
        throw(ArgumentError("--abundance must be count or rows"))

    return CliConfig(
        String(strip(input)),
        String(args["output"]),
        delimiter_char(String(args["delimiter"])),
        !args["no-default-filters"],
        filters,
        rename_count_to,
        args["min-d-region-length"],
        args["deduplicate"],
        args["use-barcode"],
        selected_threshold(args),
        linkage,
        collapse_strategy,
        soft_cutoff,
        tie_breaker,
        Float64(args["tie-atol"]),
        abundance,
        processed_output,
        args["quiet"],
    )
end

parse_filter_value(raw::AbstractString) = parse_filter_value(tryparse(Float64, raw), raw)
parse_filter_value(value::Float64, ::AbstractString) = value
parse_filter_value(::Nothing, raw::AbstractString) = raw

function make_filter(::Type{O}, column::Symbol, value::T) where {O<:AbstractFilterOperator,T}
    return ColumnFilter{O,T}(column, value)
end

function parse_filter(spec::AbstractString)
    for (operator, operator_type) in FILTER_OPERATORS
        parts = split(spec, operator; limit = 2)
        if length(parts) == 2
            column = Symbol(strip(parts[1]))
            value = parse_filter_value(strip(parts[2]))
            isempty(String(column)) && throw(ArgumentError("Invalid filter '$spec': missing column"))
            return make_filter(operator_type, column, value)
        end
    end

    throw(ArgumentError("Invalid filter '$spec'. Use operators <=, >=, ==, !=, <, or >."))
end

function resolved_column(df::AbstractDataFrame, column::Symbol)
    column in propertynames(df) && return column
    column == :V_errors && :V_error in propertynames(df) && return :V_error
    column == :J_errors && :J_error in propertynames(df) && return :J_error
    return column
end

as_number(value::Real) = Float64(value)
as_number(value) = tryparse(Float64, string(value))
as_number(::Missing) = nothing

passes_value(::LessThan, actual, expected::Real) = passes_numeric(<, actual, expected)
passes_value(::LessEqual, actual, expected::Real) = passes_numeric(<=, actual, expected)
passes_value(::GreaterThan, actual, expected::Real) = passes_numeric(>, actual, expected)
passes_value(::GreaterEqual, actual, expected::Real) = passes_numeric(>=, actual, expected)
passes_value(::EqualTo, actual, expected) = !ismissing(actual) && string(actual) == string(expected)
passes_value(::NotEqualTo, actual, expected) = !ismissing(actual) && string(actual) != string(expected)

function passes_numeric(predicate, actual, expected::Real)
    number = as_number(actual)
    return !isnothing(number) && predicate(number, Float64(expected))
end

function passes_filter(filter::ColumnFilter{O}, df::AbstractDataFrame, row::Integer) where {O<:AbstractFilterOperator}
    column = resolved_column(df, filter.column)
    return passes_value(O(), df[row, column], filter.value)
end

function apply_filter(df::DataFrame, filter::ColumnFilter; skip_missing_column::Bool)
    column = resolved_column(df, filter.column)
    if column ∉ propertynames(df)
        skip_missing_column && return df
        throw(ArgumentError("Filter column '$column' is not present in input"))
    end

    mask = [passes_filter(filter, df, row) for row in 1:nrow(df)]
    return df[mask, :]
end

function apply_filters(df::DataFrame, specs::AbstractVector{<:AbstractString}; skip_missing_column::Bool)
    result = df
    for spec in specs
        result = apply_filter(result, parse_filter(spec); skip_missing_column)
    end
    return result
end

function rename_count_column!(df::DataFrame, target::Symbol)
    :count ∈ propertynames(df) || return df
    target ∉ propertynames(df) ||
        throw(ArgumentError("Cannot rename count to '$target' because that column already exists"))
    rename!(df, :count => target)
    return df
end

function normalize_vdj_column!(df::DataFrame)
    :vdj_nt ∈ propertynames(df) && return df
    :VDJ_nt ∈ propertynames(df) && rename!(df, :VDJ_nt => :vdj_nt)
    return df
end

collapse_strategy(config::CliConfig) = collapse_strategy(Val(Symbol(config.collapse_strategy)), config)
collapse_strategy(::Val{:hardest}, ::CliConfig) = Hardest()
collapse_strategy(::Val{:soft}, config::CliConfig) = Soft(config.soft_cutoff)
collapse_strategy(::Val{:igdiscover}, ::CliConfig) = IgDiscover()

tie_breaker(config::CliConfig) = tie_breaker(Val(Symbol(replace(config.tie_breaker, "-" => "_"))))
tie_breaker(::Val{:most_common_vdj_nt}) = ByMostCommonVdjNt()
tie_breaker(::Val{:lineage_most_common_vdj_nt}) = ByLineageMostCommonVdjNt()
tie_breaker(::Val{:lowest_v_shm}) = ByLowestVShm()
tie_breaker(::Val{:vdj_count}) = ByVdjCount()
tie_breaker(::Val{:cdr3_count}) = ByCdr3Count()
tie_breaker(::Val{:sequence_count}) = BySequenceCount()
tie_breaker(::Val{:lexicographic}) = ByLexicographic()
tie_breaker(::Val{:first}) = ByFirst()
tie_breaker(::Val{:most_naive}) = ByMostNaive()

function selected_tie_breaker(config::CliConfig)
    config.collapse_strategy == "igdiscover" &&
        return ByLineageMostCommonVdjNt(; small_lineage_first_row=true)
    return tie_breaker(config)
end

function write_table(path::AbstractString, df::AbstractDataFrame; delim::Char)
    path == "-" && return CSV.write(stdout, df; delim)
    CSV.write(path, df; delim)
    return path
end

function log_progress(config::CliConfig, message::AbstractString)
    config.quiet || println(stderr, message)
    return nothing
end

function run_pipeline(config::CliConfig)
    config.quiet && return with_logger(NullLogger()) do
        run_pipeline_impl(config)
    end

    return run_pipeline_impl(config)
end

function run_pipeline_impl(config::CliConfig)
    log_progress(config, "Reading $(config.input)")
    df = CSV.read(config.input, DataFrame; delim = config.delimiter)
    log_progress(config, "Loaded $(nrow(df)) rows")

    if config.default_filters
        df = apply_filters(df, DEFAULT_FILTER_SPECS; skip_missing_column = true)
        log_progress(config, "After default filters: $(nrow(df)) rows")
    end

    df = apply_filters(df, config.filters; skip_missing_column = false)
    log_progress(config, "After custom filters: $(nrow(df)) rows")

    !isnothing(config.rename_count_to) && rename_count_column!(df, config.rename_count_to)

    preprocessed = preprocess_data(
        df;
        min_d_region_length = config.min_d_region_length,
        deduplicate = config.deduplicate,
        use_barcode = config.use_barcode,
    )
    normalize_vdj_column!(preprocessed)
    log_progress(config, "After preprocessing: $(nrow(preprocessed)) rows")

    lineages = process_lineages(preprocessed, config.threshold; linkage = config.linkage)
    nlin = nrow(lineages) == 0 ? 0 : length(unique(lineages.lineage_id))
    log_progress(config, "Processed $(nrow(lineages)) rows into $(nlin) lineages")

    if !isnothing(config.processed_output)
        write_table(config.processed_output, lineages; delim = config.delimiter)
        log_progress(config, "Wrote processed lineages to $(config.processed_output)")
    end

    collapsed = collapse_lineages(
        lineages,
        collapse_strategy(config);
        tie_breaker = selected_tie_breaker(config),
        tie_atol = config.tie_atol,
        abundance = config.abundance,
    )
    log_progress(config, "Collapsed to $(nrow(collapsed)) rows")

    write_table(config.output, collapsed; delim = config.delimiter)
    return collapsed
end

function main(argv = ARGS)
    args = parse_args(argv, build_parser())
    run_pipeline(config_from_args(args))
    return nothing
end

function julia_main()::Cint
    settings = build_parser()
    try
        main(ARGS)
        return 0
    catch err
        if err isa ArgParse.ArgParseError || err isa ArgumentError
            println(stderr, sprint(showerror, err))
            println(stderr)
            ArgParse.show_help(stderr, settings; exit_when_done=false)
            return 1
        end
        rethrow()
    end
end

end
