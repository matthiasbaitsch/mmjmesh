using GLMakie
using CairoMakie

using MMJMesh
using MMJMesh.Plots
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Utilities
using MMJMesh.Topologies
using MMJMesh.Mathematics

CairoMakie.activate!()
set_theme!(theme_minimal())
update_theme!(colormap=:jet)

meshpath(m) = joinpath(@__DIR__(), "../data/gmsh", m)
