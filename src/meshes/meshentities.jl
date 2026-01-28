"""
    MeshEntity{DT,DG,NN,G}

Entity of a finite element mesh of parametric dimension `DT` embeded in `DG` dimensional 
space connected to `NN` nodes with a geometry of type `G`.
"""
struct MeshEntity{DT,DG,NN,G}
    mesh::Mesh
    index::Integer
end

function MeshEntity(
    mesh::Mesh{DT,DG}, pdim::Integer, idx::Integer, g=GeometricObjectI{DT,DG}
) where {DT,DG}
    nn = nlinks(mesh.topology, pdim, 0, idx)
    me = MeshEntity{pdim,DG,nn,g}(mesh, idx)
    return me
end

# Short names
const Node{DG} = MeshEntity{0,DG,0}
const Edge{DG,NN} = MeshEntity{1,DG,NN}
const Face{DG,NN} = MeshEntity{2,DG,NN}
const Solid{DG,NN} = MeshEntity{3,DG,NN}

# Basic
MMJMesh.pdim(::Type{<:MeshEntity{DT}}) where {DT} = DT
MMJMesh.pdim(me::MeshEntity) = pdim(typeof(me))
MMJMesh.gdim(::Type{<:MeshEntity{DT,DG}}) where {DT,DG} = DG
MMJMesh.gdim(me::MeshEntity) = gdim(typeof(me))
MMJMesh.id(me::MeshEntity{DT}) where {DT} = MMJMesh.id(me.mesh.topology, DT, me.index)

nnodes(::Type{<:Node}) = 0
nnodes(::Type{<:MeshEntity{DT,DG,NN}}) where {DT,DG,NN} = NN
nnodes(me::MeshEntity) = nnodes(typeof(me))

geometrytype(::MeshEntity{DT,DG,NN,G}) where {DT,DG,NN,G} = G

index(me::MeshEntity) = me.index
index(me::MeshEntity{DT}, pdim::Integer, i::Integer) where {DT} =
    links(me.mesh.topology, DT, pdim)[me.index][i]
indices(me::MeshEntity{DT}, pdim::Integer) where {DT} =
    links(me.mesh.topology, DT, pdim)[me.index]
nentities(me::MeshEntity{DT}, pdim::Integer) where {DT} =
    nlinks(me.mesh.topology, DT, pdim, me.index)

entity(me::MeshEntity, pdim::Int, idx::Int) = entity(me.mesh, pdim, index(me, pdim, idx))
entities(me::MeshEntity, pdim::Integer) = MeshEntityList{pdim}(me.mesh, indices(me, pdim))

# Show
function _mpname(t::Type)
    s = string(t)
    return s[1:findfirst('{', s)-1]
end

Base.show(io::IO, e::T) where {T<:MeshEntity{0}} =
    print(io, "$(_mpname(T))($(id(e)))$(coordinates(e))")
Base.show(io::IO, e::T) where {T<:MeshEntity} =
    print(io, "$(_mpname(T))($(id(e)))$(indices(e, 0))")

# Coordinates
coordinate(n::Node, c::Int) = n.mesh.geometry.points.coordinates[c, n.index]
coordinates(n::Node) = n.mesh.geometry.points.coordinates[:, n.index]
coordinates(me::MeshEntity) = me.mesh.geometry.points.coordinates[:, indices(me, 0)]
coordinates(me::MeshEntity, i::Int) = me.mesh.geometry.points.coordinates[:, index(me, 0, i)]
coordinates(me::MeshEntity, ii::AbstractVector{Int}) =
    me.mesh.geometry.points.coordinates[:, indices(me, 0)[ii]]

# Specialized functions
Base.length(e::Edge{DG,2}) where {DG} = measure(geometry(e))
area(f::Face{2}) = measure(geometry(f))

# Geometry
geometry(node::Node{DG}) where {DG} = Point{DG}(coordinates(node))
geometry(face::Face{2,4,Box}) = Box(coordinates(face, 1), coordinates(face, 3))
geometry(me::MeshEntity{DT,DG,NN,GeometricObjectI}) where {DT,DG,NN} =
    GeometricObjectI{pdim(me)}(coordinates(me))
