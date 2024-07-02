# -------------------------------------------------------------------------------------------------
# MeshEntity struct
# -------------------------------------------------------------------------------------------------

"""
    MeshEntity{DT,DG,NN}

Entity of a finite element mesh of parametric dimension `DT` embeded in `DG` dimensional 
space connected to `NN` nodes.
"""
struct MeshEntity{DT,DG,NN}
    mesh::Mesh
    index::Int
    data::EntityData{MeshEntity{DT}}
end

function MeshEntity(mesh::Mesh{DT,DG}, pdim::Int, idx::Int) where {DT,DG}
    nn = nlinks(mesh.topology, pdim, 0, idx)
    me = MeshEntity{pdim,DG,nn}(mesh, idx, EntityData{MeshEntity{pdim}}())
    me.data.entity = me
    return me
end

# Basic
nentities(me::MeshEntity{DT}, pdim::Int) where {DT} = nlinks(me.mesh.topology, DT, pdim, me.index)
index(me::MeshEntity) = me.index
index(me::MeshEntity{DT}, pdim::Int, i::Int) where {DT} = links(me.mesh.topology, DT, pdim)[me.index][i]
indices(me::MeshEntity{DT}, pdim::Int) where {DT} = links(me.mesh.topology, DT, pdim)[me.index]
MMJMesh.pdim(::MeshEntity{DT}) where {DT} = DT
MMJMesh.gdim(::MeshEntity{DT,DG}) where {DT,DG} = DG

# Show
Base.show(io::IO, e::T) where {T<:MeshEntity} = print(io, "$(T)[$(e.index)]")

# Short names
const Node{DG} = MeshEntity{0,DG,0}
const Edge{DG,NN} = MeshEntity{1,DG,NN}
const Face{DG,NN} = MeshEntity{2,DG,NN}
const Solid{DG,NN} = MeshEntity{3,DG,NN}

# Specialized functions
Base.length(e::Edge{DG,2}) where {DG} = norm(diff(coordinates(e), dims=2))



# -------------------------------------------------------------------------------------------------
# List of entities
# -------------------------------------------------------------------------------------------------

"""
    MeshEntityList{DT}

List of entities of a finite element mesh of parametric dimension `DT`.
"""
struct MeshEntityList{DT} <: AbstractVector{MeshEntity}
    mesh::Mesh
    indices::AbstractVector{Int}
    MeshEntityList(mesh::Mesh, pdim::Int, indices::AbstractVector{Int}) = new{pdim}(mesh, indices)
end

# AbstractArray interface
Base.length(el::MeshEntityList) = length(el.indices)
Base.size(el::MeshEntityList) = (length(el),)
Base.getindex(el::MeshEntityList{DT}, i::Int) where {DT} = entity(el.mesh, DT, el.indices[i])
Base.iterate(el::MeshEntityList, state=1) = state > length(el) ? nothing : (el[state], state + 1)

# Get entities
entities(m::Mesh, pdim::Int) = MeshEntityList(m, pdim, 1:nentities(m, pdim))
entities(me::MeshEntity, pdim::Int) = MeshEntityList(me.mesh, pdim, indices(me, pdim))

# Coordinates
coordinate(n::Node, c::Int) = n.mesh.geometry.points.coordinates[c, n.index]
coordinates(n::Node) = n.mesh.geometry.points.coordinates[:, n.index]
coordinates(me::MeshEntity) = me.mesh.geometry.points.coordinates[:, indices(me, 0)]
coordinates(me::MeshEntity, i::Int) = me.mesh.geometry.points.coordinates[:, index(me, 0, i)]

# Geometry
geometry(me::MeshEntity) = GeometricObjectI{pdim(me)}(coordinates(me))

