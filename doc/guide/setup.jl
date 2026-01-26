# From Makie.jl/docs/makedocs.jl
using Pkg

cd(@__DIR__)
Pkg.activate(".")
Pkg.develop(path="../..")
Pkg.precompile()

import GLMakie
import CairoMakie

using LinearAlgebra
using CairoMakie: Figure, Axis, Axis3, scatter!, lines, lines!, DataAspect, hidespines!, hidedecorations!

using MMJMesh
using MMJMesh.Gmsh
using MMJMesh.Plots
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Utilities
using MMJMesh.Geometries
using MMJMesh.Topologies
using MMJMesh.Mathematics

CairoMakie.activate!()
CairoMakie.set_theme!(CairoMakie.theme_minimal())
CairoMakie.update_theme!(colormap=:jet)

function make_axis(f)
    ax = Axis(f, aspect=DataAspect())
    hidespines!(ax)
    hidedecorations!(ax)
    ax
end

