abstract type MeshEntity{DT, DG} end

pdim(::MeshEntity{DT, DG}) where {DT, DG} = DT
gdim(::MeshEntity{DT, DG}) where {DT, DG} = DG

"""
    nodes(me::MeshEntity)

Return the nodes of a mesh entity.
"""
nodes( ::MeshEntity) = @abstractmethod

links(me::MeshEntity{DT, DG}, d::Int) where {DT, DG} = links(me.mesh.topology, DT, d)[me.idx]
