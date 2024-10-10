using LineageCollapse
using Documenter

DocMeta.setdocmeta!(LineageCollapse, :DocTestSetup, :(using LineageCollapse); recursive=true)

makedocs(;
    modules=[LineageCollapse],
    authors="Mateusz Kaduk <mateusz.kaduk@gmail.com> and contributors",
    sitename="LineageCollapse.jl",
    format=Documenter.HTML(;
        canonical="https://mashu.github.io/LineageCollapse.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mashu/LineageCollapse.jl",
    devbranch="main",
)
