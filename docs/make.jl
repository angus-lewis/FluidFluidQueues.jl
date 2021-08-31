using Documenter, 

makedocs(
    modules = [],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "angus-lewis",
    sitename = ".jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/angus-lewis/.jl.git",
    push_preview = true
)
