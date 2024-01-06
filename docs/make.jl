push!(LOAD_PATH,"../src/")

using Documenter
using DocumenterMermaid
using MMJMesh

makedocs(
    sitename = "MMJMesh",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Home" => "index.md",
        "Meshes" => "meshes.md",
        "Associations" => "associations.md",
        "Plots" => "plots.md",
        "Utilities" => "utilities.md",
        "Geometries" => "geometries.md",
        "Topologies" => "topologies.md",
        "Gmsh" => "gmsh.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
