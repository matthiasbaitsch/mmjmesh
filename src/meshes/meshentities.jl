"""
    MeshEntity{DT,DG,NN,G}

Entity of a finite element mesh of parametric dimension `DT` embeded in `DG` dimensional 
space connected to `NN` nodes with a geometry of type `G`.
"""
struct MeshEntity{DT,DG,NN,G}
    mesh::Mesh
    index::Integer
    data::EntityData{MeshEntity{DT}}
end

function MeshEntity(
    mesh::Mesh{DT,DG}, pdim::Integer, idx::Integer, g=GeometricObjectI{DT,DG}
) where {DT,DG}
    nn = nlinks(mesh.topology, pdim, 0, idx)
    me = MeshEntity{pdim,DG,nn,g}(mesh, idx, EntityData{MeshEntity{pdim}}())
    me.data.entity = me
    return me
end

# Basic
MMJMesh.pdim(::MeshEntity{DT}) where {DT} = DT
MMJMesh.gdim(::MeshEntity{DT,DG}) where {DT,DG} = DG
nnodes(::MeshEntity{DT,DG,NN}) where {DT,DG,NN} = NN
geometrytype(::MeshEntity{DT,DG,NN,G}) where {DT,DG,NN,G} = G
nentities(me::MeshEntity{DT}, pdim::Integer) where {DT} =
    nlinks(me.mesh.topology, DT, pdim, me.index)
index(me::MeshEntity) = me.index
index(me::MeshEntity{DT}, pdim::Integer, i::Integer) where {DT} =
    links(me.mesh.topology, DT, pdim)[me.index][i]
indices(me::MeshEntity{DT}, pdim::Integer) where {DT} =
    links(me.mesh.topology, DT, pdim)[me.index]

# Show
Base.show(io::IO, e::T) where {T<:MeshEntity} = print(io, "$(T)[$(e.index)]")

# Short names
const Node{DG} = MeshEntity{0,DG,0}
const Edge{DG,NN} = MeshEntity{1,DG,NN}
const Face{DG,NN} = MeshEntity{2,DG,NN}
const Solid{DG,NN} = MeshEntity{3,DG,NN}

# Specialized functions
Base.length(e::Edge{DG,2}) where {DG} = norm(diff(coordinates(e), dims=2))

# Coordinates
coordinate(n::Node, c::Integer) = n.mesh.geometry.points.coordinates[c, n.index]
coordinates(n::Node) = n.mesh.geometry.points.coordinates[:, n.index]
coordinates(me::MeshEntity) = me.mesh.geometry.points.coordinates[:, indices(me, 0)]
coordinates(me::MeshEntity, i::Integer) = me.mesh.geometry.points.coordinates[:, index(me, 0, i)]

# Geometry
geometry(node::Node{DG}) where {DG} = Point{DG}(coordinates(node))
geometry(face::Face{2,4,Box}) = Box(coordinates(face, 1), coordinates(face, 3))
geometry(me::MeshEntity{DT,DG,NN,GeometricObjectI}) where {DT,DG,NN} =
    GeometricObjectI{pdim(me)}(coordinates(me))


"""
    MeshEntityList{DT}

List of entities of a finite element mesh of parametric dimension `DT`.
"""
struct MeshEntityList{DT} <: AbstractVector{MeshEntity}
    mesh::Mesh
    indices::AbstractVector{Int}
    MeshEntityList(mesh::Mesh, pdim::Integer, indices::AbstractVector{<:Integer}) =
        new{pdim}(mesh, indices)
end

# AbstractArray interface
Base.length(el::MeshEntityList) = length(el.indices)
Base.size(el::MeshEntityList) = (length(el),)
Base.getindex(el::MeshEntityList{DT}, i::Integer) where {DT} = entity(el.mesh, DT, el.indices[i])
Base.iterate(el::MeshEntityList, state=1) = state > length(el) ? nothing : (el[state], state + 1)

# Get entities
entities(m::Mesh, pdim::Integer) = MeshEntityList(m, pdim, 1:nentities(m, pdim)) # XXX
entities(me::MeshEntity, pdim::Integer) = MeshEntityList(me.mesh, pdim, indices(me, pdim))

