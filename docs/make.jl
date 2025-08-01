using Documenter
using Cosmos

DocMeta.setdocmeta!(Cosmos, :DocTestSetup, :(using Cosmos))

makedocs(;
	sitename = "Cosmos.jl",
	authors = "Rafael Alves Batista",
	format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
	modules = [Cosmos],
	workdir = joinpath(@__DIR__, ".."),
	pages = [
		"Home" => "index.md",
		"Installation" => "installation.md",
		"Examples" => "examples.md",
		"API Reference" => "api.md",
		],
	doctest = true,
)


deploydocs(;
	repo = "github.com/rafaelab/Cosmos.jl.git",
	branch = "gh-pages",
	versions = nothing
)