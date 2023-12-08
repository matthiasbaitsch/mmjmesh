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

# Get entity
MMJMesh.entity(m::Mesh, pdim::Int, index::Int) = MeshEntity(m, pdim, index)

# Indexes
indexes(me::MeshEntity{DT,DG,NN}, pdim::Int) where {DT,DG,NN} = links(me.mesh.topology, DT, pdim)[me.index]

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
Base.getindex(el::MeshEntityList{DT}, i::Int) where {DT} = entity(el.mesh, DT, el.indexes[i])
Base.iterate(el::MeshEntityList, state=1) = state > length(el) ? nothing : (el[state], state + 1)

# Get entities
entities(m::Mesh, pdim::Int) = MeshEntityList(m, pdim, 1:nentities(m, pdim))
entities(me::MeshEntity{DT,DG,NN}, pdim::Int) where {DT,DG,NN} = MeshEntityList(me.mesh, pdim, indexes(me, pdim))

# Coordinates
MMJMesh.coordinates(me::MeshEntity) = me.mesh.geometry.points.coordinates[:, nodeIdxs(me)]
MMJMesh.coordinates(me::MeshEntity{0, DG, NN}) where {DG,NN} = me.mesh.geometry.points.coordinates[:, me.index]
MMJMesh.coordinates(me::MeshEntity, index::Int) = me.mesh.geometry.points.coordinates[:, index]
