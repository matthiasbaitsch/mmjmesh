abstract type MeshEntity{D} end

dimension(me::MeshEntity{D}) where {D} = D

"""
    nodes(me::MeshEntity)

Return the nodes of a mesh entity.
"""
nodes(me::MeshEntity) = @abstractmethod

links(me::MeshEntity, d::Int) = links(me.mesh.topology, dimension(me), d)[me.idx]
