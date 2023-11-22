module Plots

# Exports

## plot.jl
export PlotMeshConfiguration, plot

# Modules needed by this module
import Makie
using MMJMesh.Meshes
using MMJMesh.Topologies
using MMJMesh.Geometries

# Parts
include("plot.jl")

end
