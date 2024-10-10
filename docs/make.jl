using LinageCollapse
using Documenter

DocMeta.setdocmeta!(LinageCollapse, :DocTestSetup, :(using LinageCollapse); recursive=true)

makedocs(;
    modules=[LinageCollapse],
    authors="Mateusz Kaduk <mateusz.kaduk@gmail.com> and contributors",
    sitename="LinageCollapse.jl",
    format=Documenter.HTML(;
        canonical="https://mashu.github.io/LinageCollapse.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mashu/LinageCollapse.jl",
    devbranch="main",
)
