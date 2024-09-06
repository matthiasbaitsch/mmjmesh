push!(LOAD_PATH,"../src/")

using MMJMesh
using Documenter
using DocumenterMermaid

makedocs(
    sitename = "MMJMesh",
    format = Documenter.HTML(prettyurls = false),
    build = "../build/apidoc",
    pages = [
        "Home" => "index.md",
        "Meshes" => "meshes.md",
        "Associations" => "associations.md",
        "Plots" => "plots.md",
        "Utilities" => "utilities.md",
        "Geometries" => "geometries.md",
        "Topologies" => "topologies.md",
        "Gmsh" => "gmsh.md",
        "Mathematics" => "mathematics.md"
    ]
)
