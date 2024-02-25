module Plots

# Modules needed by this module
## Others
import Makie
import MakieCore
using LinearAlgebra
using IntervalSets
## Own
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Topologies
using MMJMesh.Mathematics

# Exports
export mplot, mconf

# Parts
include("sampleadaptive.jl")
include("mplot.jl")
include("fplot.jl")

end
