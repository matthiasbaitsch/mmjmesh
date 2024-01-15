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

    MeshEntity(mesh::Mesh{DT,DG}, pdim::Int, index::Int) where {DT,DG} =
        new{pdim,DG,length(links(mesh.topology, pdim, 0), index)}(mesh, index)
end

# Names
const Node{DG} = MeshEntity{0,DG,0}
const Edge{DG,NN} = MeshEntity{1,DG,NN}
const Face{DG,NN} = MeshEntity{2,DG,NN}
const Solid{DG,NN} = MeshEntity{2,DG,NN}

# Basic
entity(m::Mesh, pdim::Int, index::Int) = MeshEntity(m, pdim, index)
nentities(me::MeshEntity{DT,DG,NN}, pdim::Int) where {DT,DG,NN} = nlinks(me.mesh.topology, DT, pdim, me.index)
index(me::MeshEntity{DT,DG,NN}, pdim::Int, i::Int) where {DT,DG,NN} = links(me.mesh.topology, DT, pdim)[me.index][i]
indexes(me::MeshEntity{DT,DG,NN}, pdim::Int) where {DT,DG,NN} = links(me.mesh.topology, DT, pdim)[me.index]
Base.length(e::Edge{DG,2}) where {DG} = norm(diff(coordinates(e), dims=2))

# Show
name(e::MeshEntity) = name(typeof(e))
name(::Type{MeshEntity{0,DG,NN}}) where {DG,NN} = "Node"
name(::Type{MeshEntity{1,DG,NN}}) where {DG,NN} = "Edge"
name(::Type{MeshEntity{2,DG,NN}}) where {DG,NN} = "Face"
name(::Type{MeshEntity{3,DG,NN}}) where {DG,NN} = "Solid"
Base.show(io::IO, e::MeshEntity) = print(io, "$(name(e))[$(e.index)]")


# -------------------------------------------------------------------------------------------------
# List of entities
# -------------------------------------------------------------------------------------------------

"""
    MeshEntityList{DT}

List of entities of a finite element mesh of parametric dimension `DT`.
"""
struct MeshEntityList{DT} <: AbstractVector{MeshEntity}
    mesh::Mesh
    indexes::AbstractVector{Int}
    MeshEntityList(mesh::Mesh, pdim::Int, indexes::AbstractVector{Int}) = new{pdim}(mesh, indexes)
end

# AbstractArray interface
Base.length(el::MeshEntityList) = length(el.indexes)
Base.size(el::MeshEntityList) = (length(el),)
Base.getindex(el::MeshEntityList{DT}, i::Int) where {DT} = entity(el.mesh, DT, el.indexes[i])
Base.iterate(el::MeshEntityList, state=1) = state > length(el) ? nothing : (el[state], state + 1)

# Get entities
entities(m::Mesh, pdim::Int) = MeshEntityList(m, pdim, 1:nentities(m, pdim))
entities(me::MeshEntity{DT,DG,NN}, pdim::Int) where {DT,DG,NN} = MeshEntityList(me.mesh, pdim, indexes(me, pdim))

# Coordinates
coordinate(n::Node, c::Int) = n.mesh.geometry.points.coordinates[c, n.index]
coordinates(n::Node) = n.mesh.geometry.points.coordinates[:, n.index]

coordinates(me::MeshEntity) = me.mesh.geometry.points.coordinates[:, indexes(me, 0)]
coordinates(me::MeshEntity, i::Int) = me.mesh.geometry.points.coordinates[:, index(me, 0, i)]
