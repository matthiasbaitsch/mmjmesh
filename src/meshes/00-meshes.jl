module Meshes

# Exports

## mesh.jl
export Mesh
export pdim, gdim

## meshentities.jl
export Node, Edge, Face, Solid
export dimension, nodes, nodeindexes

# Modules needed by this module
using MMJMesh.MMJBase
using MMJMesh.Topologies

# Parts
include("mesh.jl")
include("meshentities.jl")

end