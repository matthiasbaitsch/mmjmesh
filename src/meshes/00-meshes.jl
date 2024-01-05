module Meshes

# Modules needed by this module
## Other
using LinearAlgebra
import Base.length
import Base.size
## Own
using MMJMesh
using MMJMesh.MMJBase
using MMJMesh.Topologies
using MMJMesh.Geometries
using MMJMesh.Associations

# Exports
## mesh.jl
export Mesh
export pdim, gdim, nelements, element, elements
## meshentities.jl
export MeshEntity
## common.jl
export coordinates
export nnodes, nedges, nfaces, nsolids
export node, edge, face, solid
export nodes, edges, faces, solids
export nodeIdxs, edgeIdxs, faceIdxs, solidIdxs

# Parts
include("mesh.jl")
include("meshentities.jl")
include("common.jl")

end