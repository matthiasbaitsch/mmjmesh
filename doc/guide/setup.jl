import GLMakie
import CairoMakie
import CairoMakie: Figure, Axis3, scatter!, lines

using MMJMesh
using MMJMesh.Plots
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Utilities
using MMJMesh.Topologies
using MMJMesh.Mathematics

using DomainSets: ×, (..)

CairoMakie.activate!()
CairoMakie.set_theme!(CairoMakie.theme_minimal())
CairoMakie.update_theme!(colormap=:jet)

meshpath(m) = joinpath(@__DIR__(), "../../data/gmsh", m)
