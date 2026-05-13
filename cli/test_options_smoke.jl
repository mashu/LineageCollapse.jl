#!/usr/bin/env julia
# Smoke-test CLI flag combinations. Default INPUT is a small slice of real data for speed.
# Example: julia --project=cli cli/test_options_smoke.jl /path/to/file.tsv.gz

using LineageCollapseCLI

function smoke_main(argv::Vector{String})
    isempty(argv) && error("usage: test_options_smoke.jl INPUT.tsv.gz")
    input = argv[1]
    isfile(input) || error("input not found: $input")

    failures = Tuple{String,String}[]
    function run_case(label::AbstractString, args::Vector{String})
        try
            LineageCollapseCLI.main(args)
            println("OK  ", label)
        catch err
            msg = sprint(showerror, err)
            push!(failures, (label, msg))
            println("FAIL ", label, ": ", msg)
        end
    end

    id = Ref(0)
    next_id() = (id[] += 1; id[])
    out() = joinpath(tempdir(), "lc-smoke-$(next_id()).tsv")
    proc_path() = joinpath(tempdir(), "lc-smoke-proc-$(next_id()).tsv")

    base = [input, "--quiet"]

    run_case("defaults", vcat(base, ["-o", out()]))
    run_case("no-default-filters", vcat(base, ["--no-default-filters", "-o", out()]))
    run_case("filter productive", vcat(base, ["--filter", "productive==true", "-o", out()]))
    run_case("filter junction_length", vcat(base, ["--filter", "junction_length>=10", "-o", out()]))
    run_case("keep-count", vcat(base, ["--keep-count", "-o", out()]))
    run_case("rename-count-to", vcat(base, ["--rename-count-to", "my_reads", "-o", out()]))
    run_case("min-d-region-length 0", vcat(base, ["--min-d-region-length", "0", "-o", out()]))
    run_case("deduplicate", vcat(base, ["--deduplicate", "-o", out()]))
    run_case("deduplicate+use-barcode", vcat(base, ["--deduplicate", "--use-barcode", "-o", out()]))
    run_case("cdr3-differences 0", vcat(base, ["--cdr3-differences", "0", "-o", out()]))
    run_case("cdr3-differences 2", vcat(base, ["--cdr3-differences", "2", "-o", out()]))
    run_case("cdr3-distance-fraction", vcat(base, ["--cdr3-distance-fraction", "0.12", "-o", out()]))
    run_case("cdr3-similarity-percent", vcat(base, ["--cdr3-similarity-percent", "88", "-o", out()]))
    run_case("linkage complete", vcat(base, ["--linkage", "complete", "-o", out()]))
    run_case("linkage average", vcat(base, ["--linkage", "average", "-o", out()]))
    run_case("collapse soft", vcat(base, ["--collapse", "soft", "--soft-cutoff", "0.4", "-o", out()]))
    pout = proc_path()
    run_case("processed-output", vcat(base, ["--processed-output", pout, "-o", out()]))

    for tb in (
        "most-common-vdj-nt",
        "vdj-count",
        "cdr3-count",
        "sequence-count",
        "lexicographic",
        "first",
        "most-naive",
    )
        run_case("tie-breaker $tb", vcat(base, ["--tie-breaker", tb, "-o", out()]))
    end

    run_case("tie-atol", vcat(base, ["--tie-atol", "0.02", "-o", out()]))

    # delimiter: IgBLAST is tab; explicit "\\t" in ArgParse argv is the string \t one char when parsed by shell - here we pass tab via Char
    run_case("delimiter tab literal", vcat(base, ["--delimiter", "\t", "-o", out()]))

    # Conflicting thresholds should error before heavy work
    try
        LineageCollapseCLI.main(vcat(base, [
            "--cdr3-distance-fraction", "0.1",
            "--cdr3-similarity-percent", "90",
            "-o", out(),
        ]))
        push!(failures, ("conflict distance+similarity", "expected error"))
        println("FAIL conflict distance+similarity: expected error")
    catch err
        println("OK  conflict distance+similarity throws")
    end

    if !isempty(failures)
        println(stderr, "\n$(length(failures)) failure(s):")
        for (label, msg) in failures
            println(stderr, "  - ", label, ": ", msg)
        end
        return 1
    end
    return 0
end

exit(smoke_main(ARGS))
