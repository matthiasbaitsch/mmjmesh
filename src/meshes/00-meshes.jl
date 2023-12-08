module Meshes

# Modules needed by this module
using MMJMesh
using MMJMesh.MMJBase
using MMJMesh.Topologies

"""
    nnodes(m)

Number of nodes associated with `m`.
"""
nnodes(m) = nentities(m, 0)

"""
    nedges(m)

Number of edges associated with `m`.
"""
nedges(m) = nentities(m, 1)

"""
    nfaces(m)

Number of faces associated with `m`.
"""
nnfaces(m) = nentities(m, 2)

"""
    nsolids(m)

Number of solids associated with `m`.
"""
nsolids(m) = nentities(m, 3)

"""
    nodeIdxs(m)

Indexes of nodes associated with `m`.
"""
nodeIdxs(m) = indexes(m, 0)

"""
    edgeIdxs(m)

Indexes of edges associated with `m`.
"""
edgeIdxs(m) = indexes(m, 1)

"""
    faceIdxs(m)

Indexes of faces associated with `m`.
"""
faceIdxs(m) = indexes(m, 2)

"""
    solidIdxs(m)

Indexes of solids associated with `m`.
"""
solidIdxs(m) = indexes(m, 3)

"""
    node(m, index::Int)

Node of `m` at `index`.
"""
node(m, index::Int) = entity(m, 0, index)

"""
    edge(m, index::Int)

Edge of `m` at `index`.
"""
edge(m, index::Int) = entity(m, 1, index)

"""
    face(m, index::Int)

Face of `m` at `index`.
"""
face(m, index::Int) = entity(m, 2, index)

"""
    solid(m, index::Int)

Solid of `m` at `index`.
"""
solid(m, index::Int) = entity(m, 3, index)

"""
    nodes(m)

Nodes of `m`.
"""
nodes(m) = entities(m, 0)

"""
    edges(m)

Edges of `m`.
"""
edges(m) = entities(m, 1)

"""
    faces(m)

Faces of `m`.
"""
faces(m) = entities(m, 2)

"""
    solids(m)

Solids of `m`.
"""
solids(m) = entities(m, 3)


# Exports
## mesh.jl
export Mesh
export pdim, gdim
## meshentities.jl
export MeshEntity
export node, edge, face, solid
export nodes, edges, faces, solids
export nodeIdxs, edgeIdxs, faceIdxs, solidIdxs

# Parts
include("mesh.jl")
include("meshentities.jl")

end