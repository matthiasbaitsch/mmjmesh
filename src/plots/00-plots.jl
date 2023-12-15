module Plots

# Modules needed by this module
## Others
import Makie
using LinearAlgebra
## Own
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Topologies

# Exports
export mplot, mconf

# Parts
include("plotrecipe.jl")
include("styles.jl")
include("doplot.jl")

end
