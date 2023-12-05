module Plots

# Exports

## plot.jl
export PlotStyle, plot

# Modules needed by this module
import Makie
using LinearAlgebra
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Topologies
using MMJMesh.Geometries

# Parts
include("style.jl")
include("plot.jl")

end
