# -------------------------------------------------------------------------------------------------
# MeshEntity struct
# -------------------------------------------------------------------------------------------------

"""
    Mesh(dt, dg [, nn=0])
    Mesh(coordinates, dt)
    Mesh(coordinates, elements, dt)

Mesh of parametric dimension ``d_t`` embedded in ``d_g`` dimensional space.
"""
struct Mesh{DT,DG}
    topology::Topology{DT}
    geometry::Geometry{DG}
    groups::EntityGroupCollection
    data::MeshData
end

Mesh(dt::Int, dg::Int, nn::Int=0) = Mesh{dt,dg}(Topology(dt, nn), Geometry(dg, nn), EntityGroupCollection(), MeshData(Mesh{dt,dg}))

Mesh(coordinates::Matrix, dt::Int) =
    Mesh{dt,size(coordinates, 1)}(
        Topology(dt, size(coordinates, 2)),
        Geometry(coordinates),
        EntityGroupCollection(),
        MeshData(Mesh{dt,size(coordinates, 1)})
    )

function Mesh(coordinates::Matrix, elements::Vector{Vector{Int}}, dt::Int)
    m = Mesh(coordinates, dt)
    addlinks!(m.topology, dt, 0, elements)
    populatepredfinedgroups!(m)
    return m
end

nentities(m::Mesh, dim::Int) = Topologies.nentities(m.topology, dim, true)
indexes(m::Mesh, pdim::Int) = 1:nentities(m, pdim)
pdim(::Mesh{DT, DG}) where {DT, DG} = DT
gdim(::Mesh{DT, DG}) where {DT, DG} = DG
coordinates(m::Mesh) = m.geometry.points.coordinates[:, nodeIdxs(m)]
coordinates(m::Mesh, index::Int) = m.geometry.points.coordinates[:, index]

nelements(m::Mesh{DT, DG}) where {DT, DG} = nentities(m, DT)
element(m::Mesh{DT, DG}, index::Int) where {DT, DG} = entity(m, DT, index)
elements(m::Mesh{DT, DG}) where {DT, DG} = entities(m, DT)

function Base.show(io::IO, m::Mesh{DT,DG}) where {DT,DG}
    n = 100
    s1 = repeat("=", n)
    s2 = repeat("-", n)
    println(io, s1)
    println(io, "Mesh{", DT, ", ", DG, "}")
    println(io, s1)
    print(io, m.topology)
    println(io, s2)
    print(io, m.geometry)
    println(io, s2)
    print(io, m.groups)
    println(io, s1)
end


