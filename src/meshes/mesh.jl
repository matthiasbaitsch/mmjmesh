using MMJMesh.Topologies: Topology
using MMJMesh.Geometries: Geometry

struct Mesh{DT,DG}
    topology::Topology{DT}
    geometry::Geometry{DG}
end

Mesh(dt::Int, dg::Int, nn::Int=0) = Mesh{dt,dg}(Topology(dt, nn), Geometry(dg, nn))

Mesh(coordinates::Matrix{T}, dt::Int=size(coordinates, 1)) where {T<:Number} =
    Mesh{dt,size(coordinates, 1)}(
        Topology(dt, size(coordinates, 2)),
        Geometry(coordinates))

nentities(m::Mesh, dim::Int) = Topologies.nentities(m.topology, dim, true)
indexes(m::Mesh, pdim::Int) = 1:nentities(m, pdim)
pdim(m::Mesh) = Topologies.dimension(m.topology)
gdim(m::Mesh) = Geometries.dimension(m.geometry)
MMJMesh.coordinates(m::Mesh) = m.geometry.points.coordinates[:, nodeIdxs(m)]
MMJMesh.coordinates(m::Mesh, index::Int) = m.geometry.points.coordinates[:, index]

function Base.show(io::IO, m::Mesh{DT,DG}) where {DT,DG}
    n = 100
    s1 = repeat("=", n)
    s2 = repeat("-", n)
    println(io, s1)
    println(io, "Mesh{", DT, ", ", DG, "}")
    println(io, s1)
    show(io, m.topology)
    println(io, s2)
    show(io, m.geometry)
    println(io, s1)
end


