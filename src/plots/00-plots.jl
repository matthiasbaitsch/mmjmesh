module Plots

# Modules needed by this module

## Others
using Random
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
include("sampleadaptive2.jl")
include("mplot.jl")
include("fplot.jl")

end
