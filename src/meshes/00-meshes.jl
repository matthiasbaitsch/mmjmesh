module Meshes

# Modules needed by this module
using MMJMesh
using MMJMesh.MMJBase
using MMJMesh.Topologies

# Exports
## mesh.jl
export Mesh
export pdim, gdim
## meshentities.jl
export MeshEntity
## common.jl
export nnodes, nedges, nfaces, nsolids
export node, edge, face, solid
export nodes, edges, faces, solids
export nodeIdxs, edgeIdxs, faceIdxs, solidIdxs

# Parts
include("mesh.jl")
include("meshentities.jl")
include("common.jl")

end