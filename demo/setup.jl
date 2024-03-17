import CairoMakie as cm

using MMJMesh
using MMJMesh.Plots
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Utilities
using MMJMesh.Topologies
using MMJMesh.Mathematics

cm.set_theme!(cm.theme_minimal())
cm.update_theme!(colormap=:jet)

meshpath(m) = joinpath(@__DIR__(), "../data/gmsh", m)
