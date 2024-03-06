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
using MMJMesh.Topologies
using MMJMesh.Mathematics

# Exports
export mplot, mconf

# Parts
include("curveapproximation.jl")
include("sampleadaptive.jl")
include("mplot.jl")
include("fplot.jl")

end
