"""
    MeshEntityList{DT}

List of entities of a finite element mesh of parametric dimension `DT`.
"""
struct MeshEntityList{DT} <: AbstractVector{MeshEntity{DT}}
    mesh::Mesh
    indices::AbstractVector{<:Integer}
end

# Misc
mesh(mel::MeshEntityList) = mel.mesh
edim(::MeshEntityList{DT}) where DT = DT

# AbstractArray interface
Base.size(el::MeshEntityList) = (length(el),)
Base.length(el::MeshEntityList) = length(el.indices)
Base.getindex(el::MeshEntityList{DT}, i::Integer) where DT = entity(el.mesh, DT, el.indices[i])
Base.iterate(el::MeshEntityList, state=1) = state > length(el) ? nothing : (el[state], state + 1)

# Indices
indices(mel::MeshEntityList) = mel.indices

# Get entities from mesh entity list
entities(mel::MeshEntityList) = mel
entities(mel::MeshEntityList, pdim::Integer) = entities(mesh(mel), pdim, indices(mel, pdim))

