module Plots

# Exports

## plot.jl
export PlotMeshConfiguration, plotMesh

# Modules needed by this module
using MMJMesh.Meshes
using MMJMesh.Topologies
using MMJMesh.Geometries
# using Makie: scatter!, hidedecorations!, Figure, Axis, DataAspect

# Parts
include("plot.jl")

end
