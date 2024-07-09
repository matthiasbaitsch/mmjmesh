module Plots

# Modules needed by this module

## Others
using Printf
using Random
using StaticArrays
using IntervalSets
using LinearAlgebra
using DomainSets: Rectangle

import Makie
import MakieCore

## Own
using MMJMesh
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Utilities
using MMJMesh.Geometries
using MMJMesh.Topologies
using MMJMesh.Mathematics

# Exports
export mplot, mconf, fplot, fplot3d, vplot

# TODO good place or better, get rid of
function _getcolor(x::Matrix, color, zscale)
    if typeof(color) == Int && 1 <= color <= 3
        return x[color, :] / zscale
    end
    return color
end

# Parts
include("helperfunctions.jl")
include("curveapproximation.jl")
include("sample1d.jl")
include("sample2d.jl")
include("mplot.jl")
include("fplot.jl")
include("vplot.jl")
include("fplot3d.jl")

end
