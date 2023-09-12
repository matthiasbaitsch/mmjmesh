using MMJMesh.Topologies: Topology
using MMJMesh.Geometries: Geometry

import Base.show

struct Mesh{DT,DG}
    topology::Topology{DT}
    geometry::Geometry{DG}
end

Mesh(dt::Int, dg::Int, nn::Int=0) = Mesh{dt,dg}(Topology(dt, nn), Geometry(dg, nn))

Mesh(coordinates::Matrix{T}, dt::Int=size(coordinates, 1)) where T<:Number =
    Mesh{dt,size(coordinates, 1)}(
        Topology(dt, size(coordinates, 2)),
        Geometry(coordinates))

function show(io::IO, m::Mesh{DT,DG}) where {DT,DG}
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


