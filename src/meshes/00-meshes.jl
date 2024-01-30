module Meshes

# Modules needed by this module
using LinearAlgebra

using MMJMesh
using MMJMesh.MMJBase
using MMJMesh.Topologies
using MMJMesh.Geometries

# Exports
## groups.jl
export Group, GroupCollection, addrecipe!, clearcache!, names
## mesh.jl
export Mesh
export pdim, gdim, nelements, element, elements, entity, entities
## meshentities.jl
export MeshEntity, Node, Edge, Face, Solid, index
## common.jl
export coordinate, coordinates
export nnodes, nedges, nfaces, nsolids
export node, edge, face, solid
export nodes, edges, faces, solids
export nodeindices, edgeindices, faceindices, solidindices
## entitygroups.jl
export collectgroups, groupids, populatepredfinedgroups!, edim
export ngroups, groupnames, hasgroups
export groupof, groupsof
export NodeGroup, EdgeGroup, FaceGroup, SolidGroup

# Parts
include("data.jl")
include("groups.jl")
include("mesh.jl")
include("entitydata.jl")
include("meshentities.jl")
include("common.jl")
include("entitygroups.jl")

end