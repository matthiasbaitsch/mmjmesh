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

## Own
using MMJMesh
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Geometries
using MMJMesh.Topologies
using MMJMesh.Mathematics

# Exports
export mplot, fplot, fplot3d, mconf

# Parts
include("curveapproximation.jl")
include("sample1d.jl")
include("sample2d.jl")
include("mplot.jl")
include("fplot.jl")
include("fplot3d.jl")

# TODO
# plottype(::T) where {T<:AbstractMapping{SVector{2,Float64},Real}} = FPlot3D

end
