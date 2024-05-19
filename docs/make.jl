using Documenter

push!(LOAD_PATH, "..")

using Cosmos

DocMeta.setdocmeta!(Cosmos, :DocTestSetup, :(using Cosmos))


makedocs(
	sitename = "Cosmos.jl",
	format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
	modules = [Cosmos],
	workdir = joinpath(@__DIR__, ".."),
	pages = [
		"Home" => "index.md",
	]
)

deploydocs(repo = "github.com/rafaelab/Cosmos.jl.git", versions = nothing)