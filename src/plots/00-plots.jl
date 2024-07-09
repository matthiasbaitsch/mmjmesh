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
export mplot, mconf, fplot, fplot3d, vplot

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
