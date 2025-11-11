# From Makie.jl/docs/makedocs.jl
using Pkg

cd(@__DIR__)
Pkg.activate(".")
Pkg.develop(path="../..")
Pkg.precompile()

import GLMakie
import CairoMakie

using LinearAlgebra
using CairoMakie: Figure, Axis3, scatter!, lines, lines!

using MMJMesh
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

meshpath(m) = joinpath(@__DIR__(), "../../data/gmsh", m);
