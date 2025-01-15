"""
    Mesh(dt, dg [, nn=0, g1=GeometricObjectI, g2=GeometricObjectI])
    Mesh(coordinates, dt [, g1=GeometricObjectI, g2=GeometricObjectI])
    Mesh(coordinates, elements, dt [, g1=GeometricObjectI, g2=GeometricObjectI])

Mesh of parametric dimension `dt` embedded in `dg` dimensional space. Types `g1` and `g2` 
are default types for geometric objects of parametric dimension 1 and 2.
"""
struct Mesh{DT,DG,G1,G2}
    topology::Topology{DT}
    geometry::Geometry{DG}
    groups::GroupCollection
    data::Data
end

function Mesh(
    dt::Int, dg::Int, nn::Int=0;
    g1=GeometricObjectI, g2=GeometricObjectI
)
    T = Mesh{dt,dg,g1,g2}
    mesh = T(Topology(dt, nn), Geometry(dg, nn), GroupCollection(), Data())
    mesh.data.mesh = mesh
    _registergrouprecipies(mesh)
    return mesh
end

function Mesh(
    coordinates::Matrix, dt::Int;
    g1=GeometricObjectI, g2=GeometricObjectI
)
    m = Mesh(dt, size(coordinates, 1), size(coordinates, 2), g1=g1, g2=g2)
    m.geometry.points.coordinates[:, :] = coordinates[:, :] # TODO think again
    return m
end

function Mesh(
    coordinates::Matrix,
    elements::Vector{Vector{Int}}, dt::Int;
    g1=GeometricObjectI, g2=GeometricObjectI
)
    m = Mesh(coordinates, dt, g1=g1, g2=g2)
    addlinks!(m.topology, dt, 0, elements)
    return m
end

nentities(m::Mesh, dim::Int) = Topologies.nentities(m.topology, dim, true)
entity(m::Mesh, pdim::Int, idx::Int) = MeshEntity(m, pdim, idx, geometrytype(m, pdim))
indices(m::Mesh, pdim::Int) = 1:nentities(m, pdim)

# Geometry
function geometrytype(::Mesh{DT,DG,G1,G2}, pdim::Integer) where {DT,DG,G1,G2}
    pdim == 0 && return Point
    pdim == 1 && return G1
    pdim == 2 && return G2
    error()
end

MMJMesh.pdim(::Mesh{DT}) where {DT} = DT
MMJMesh.gdim(::Mesh{DT,DG}) where {DT,DG} = DG

coordinates(m::Mesh) = m.geometry.points.coordinates[:, nodeindices(m)]
coordinates(m::Mesh, group::Symbol) = m.geometry.points.coordinates[:, m.groups[group]]
coordinates(m::Mesh, index::Int) = m.geometry.points.coordinates[:, index]

nelements(m::Mesh{DT,DG}) where {DT,DG} = nentities(m, DT)
element(m::Mesh{DT,DG}, index::Int) where {DT,DG} = entity(m, DT, index)
elements(m::Mesh{DT,DG}) where {DT,DG} = entities(m, DT)

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


