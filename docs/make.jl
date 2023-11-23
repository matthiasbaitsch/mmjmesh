push!(LOAD_PATH,"../src/")

using Documenter
using MMJMesh

makedocs(
    sitename = "MMJMesh",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Home" => "index.md",
        "Topologies" => "topologies.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
