#!/usr/bin/env julia

using Pkg

const CLI_DIR = @__DIR__
const PACKAGE_ROOT = dirname(CLI_DIR)
const DEFAULT_APP_DIR = joinpath(PACKAGE_ROOT, "build", "lineage-collapse")

function argument_value(flag::String, default::String)
    index = findfirst(==(flag), ARGS)
    isnothing(index) && return default
    index < length(ARGS) || throw(ArgumentError("Missing value for $flag"))
    return ARGS[index + 1]
end

function flag_present(flag::String)
    return flag in ARGS
end

function prepare_environment()
    Pkg.activate(CLI_DIR)
    Pkg.develop(PackageSpec(path = PACKAGE_ROOT))
    Pkg.instantiate()
    return nothing
end

prepare_environment()

using PackageCompiler

function build_app()
    app_dir = argument_value("--output", DEFAULT_APP_DIR)
    cpu_target = argument_value("--cpu-target", "generic")
    force = flag_present("--force")

    create_app(
        CLI_DIR,
        app_dir;
        executables = ["lineage-collapse" => "julia_main"],
        precompile_execution_file = joinpath(CLI_DIR, "precompile.jl"),
        force,
        cpu_target,
        filter_stdlibs = true,
        include_lazy_artifacts = false,
        sysimage_build_args = `-O3`,
    )

    println("Built lineage-collapse app at $app_dir")
    println("Executable: $(joinpath(app_dir, "bin", "lineage-collapse"))")
    return nothing
end

build_app()
