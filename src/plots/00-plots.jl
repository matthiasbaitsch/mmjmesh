module Plots

# Modules needed by this module
## Others
import Makie
using LinearAlgebra
## Own
using MMJMesh
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Topologies
using MMJMesh.Geometries

# Exports
## styles.jl
export PlotStyle
## plot.jl
export plot

# Parts
include("styles.jl")
include("plot.jl")

end
