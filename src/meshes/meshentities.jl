abstract type MeshEntity{D} end

dimension(me::MeshEntity{D}) where {D} = D

"""
    nodes(me::MeshEntity)

Return the nodes of a mesh entity.
"""
nodes(me::MeshEntity) = @abstractmethod

export links
links(me::MeshEntity, d::Int) = links(me.mesh.topology, dimension(me), d)[me.idx]

struct Node <: MeshEntity{0}
    mesh::Mesh
    idx::Int
end

struct Edge <: MeshEntity{1}
    mesh::Mesh
    idx::Int
end

struct Face <: MeshEntity{2}
    mesh::Mesh
    idx::Int
end

struct Solid <: MeshEntity{3}
    mesh::Mesh
    idx::Int
end

