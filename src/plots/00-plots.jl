module Plots

# Modules needed by this module

## Others
using Printf
using StaticArrays
using LinearAlgebra

import Makie
import Random
import MakieCore
import DomainSets
import IntervalSets

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

end
