using Pkg
Pkg.add("Documenter")
using Documenter
Pkg.activate(".")
# push!(LOAD_PATH,"./src/")

using Documenter, dcporbit

makedocs(sitename="Discontinuous particle orbit solver",
    pages = [
        "Home" => "index.md"
    ],
    format=Documenter.HTML(prettyurls=false),
    modules = [dcporbit]
    )

deploydocs(
    repo = "github.com/Spiffmeister/dcporbit.jl.git"
)