module Plots

# Modules needed by this module
import Makie
using LinearAlgebra
using MMJMesh
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Topologies
using MMJMesh.Geometries

# Exports
## plot.jl
export PlotStyle, plot

# Parts
include("style.jl")
include("plot.jl")

end
