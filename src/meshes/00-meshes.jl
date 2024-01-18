module Meshes

# Modules needed by this module
using LinearAlgebra

using MMJMesh
using MMJMesh.MMJBase
using MMJMesh.Topologies
using MMJMesh.Geometries

# Exports
## groups.jl
export EntityGroup, NodeGroup, EdgeGroup, FaceGroup, SolidGroup
export ngroups, groupnames, hasgroups
## mesh.jl
export Mesh
export pdim, gdim, nelements, element, elements, entity, entities
## meshentities.jl
export MeshEntity, Node, Edge, Face, Solid
## common.jl
export coordinate, coordinates
export nnodes, nedges, nfaces, nsolids
export node, edge, face, solid
export nodes, edges, faces, solids
export nodeIdxs, edgeIdxs, faceIdxs, solidIdxs
## entitygroups.jl
export collectgroups, groupids, populatepredfinedgroups!, edim

# Parts
include("data.jl")
include("groups.jl")
include("mesh.jl")
include("entitydata.jl")
include("meshentities.jl")
include("common.jl")
include("entitygroups.jl")

end