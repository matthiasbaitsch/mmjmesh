module Utilities

# Exports
## generatemeshes.jl
export makemeshonrectangle
## plot.jl
export PlotMeshConfiguration, plotMesh

# Module needed by this module
using MMJMesh.Meshes
using MMJMesh.Topologies
using MMJMesh.Geometries

using Makie: scatter!, hidedecorations!, Figure, Axis, DataAspect

# Parts
include("generatemeshes.jl")
include("plot.jl")

end
