module Meshes

export Mesh
export Node, Edge, Face, Solid

export dimension, nodes, nodeindexes

using MMJMesh.MMJBase
using MMJMesh.Topologies

include("mesh.jl")

include("meshentities.jl")

end