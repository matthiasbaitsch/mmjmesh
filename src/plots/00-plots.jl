module Plots

# Modules needed by this module

## Others
using Printf
using Random
using StaticArrays
using IntervalSets
using LinearAlgebra

import Makie
import MakieCore

using DomainSets: Rectangle, component, components

## Own
using MMJMesh
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Utilities
using MMJMesh.Geometries
using MMJMesh.Topologies
using MMJMesh.Mathematics

# Exports
export plot, plot!, mplot, mconf, fplot, fplot3d, vplot, feplot, feconf

# Parts
include("helperfunctions.jl")
include("curveapproximation.jl")
include("sample1d.jl")
include("sample2d.jl")
include("mplot.jl")
include("fplot.jl")
include("vplot.jl")
include("fplot3d.jl")
include("feplot.jl")
include("symbols.jl")

plot(m::Mesh, args...; kwargs...) = mplot(m, args...; kwargs...) |> mconf()
plot!(m::Mesh, args...; kwargs...) = mplot!(m, args...; kwargs...) |> mconf()
plot(f::AbstractMapping, args...; kwargs...) = fplot(f, args...; kwargs...)
plot!(f::AbstractMapping, args...; kwargs...) = fplot!(f, args...; kwargs...)

end
