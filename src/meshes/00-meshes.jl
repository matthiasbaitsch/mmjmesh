module Meshes

# Modules needed by this module
using LinearAlgebra

using MMJMesh
using MMJMesh.MMJBase
using MMJMesh.Topologies
using MMJMesh.Geometries

import MMJMesh.MMJBase: gdim, pdim

# Exports
## groups.jl
export Group, GroupCollection, addrecipe!, clearcache!, names
## mesh.jl
export Mesh
export nelements, element, elements, entity, entities
## meshentities.jl
export MeshEntity, Node, Edge, Face, Solid, index, geometry
## common.jl
export coordinate, coordinates
export nnodes, nedges, nfaces, nsolids
export node, edge, face, solid
export nodes, edges, faces, solids
export nodeindices, edgeindices, faceindices, solidindices
## entitygroups.jl
export collectgroups, groupids, populatepredfinedgroups!, edim
export ngroups, groupnames, hasgroups
export group, groups
export NodeGroup, EdgeGroup, FaceGroup, SolidGroup

# Parts
include("data.jl")
include("groups.jl")
include("mesh.jl")
include("meshentities.jl")
include("entitygroups.jl")
include("entitydata.jl")
include("common.jl")

end