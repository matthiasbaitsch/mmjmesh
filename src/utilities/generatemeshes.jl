_isscaling(gmap) = (gmap == identity)
_isscaling(m::AffineMapping) = m.A[1, 2] == 0 && m.A[2, 1] == 0


function makemeshoninterval(a::Number, b::Number, n::Int, g=t -> [t; 0.0])
    # Coordinates
    c = reduce(hcat, g.(range(a, stop=b, length=n + 1)))

    # Mesh
    m = Mesh(c, 1)

    # Coordinates
    nn = n + 1
    c = zeros(2, nn)
    c[1, :] = range(start=1, stop=b, length=nn)

    # Connectivity
    cl = ConnectivityList()
    for i in 1:n
        push!(cl, [i, i + 1])
    end
    addlinks!(m.topology, 1, 0, cl)

    # Done
    return m
end

"""
    Mesh(I::Interval, n::Int)

Construct uniform mesh by dividing interval `I` into `n` elements.
"""
Meshes.Mesh(I::IntervalSets.AbstractInterval, n::Int, g=t -> [t; 0.0]) =
    makemeshoninterval(IntervalSets.leftendpoint(I), IntervalSets.rightendpoint(I), n, g)

"""
    Meshtype

Types of meshes on structured grids.
"""
@enum Meshtype begin
    QUADRANGLE
    TRIANGLE
end

"""
Generate a quad mesh on the specified with `nex` elements in the 
``x`` direction and `ney` elements in the ``y`` direction.
"""
makemeshonrectangle(
    w::Real, h::Real, nex::Integer, ney::Integer=-1, meshtype::Meshtype=QUADRANGLE; gmap=identity
) = makemeshonrectangle((0 .. w) × (0 .. h), nex, ney, meshtype, gmap=gmap)

function makemeshonrectangle(
    Ω::DomainSets.Rectangle, nex::Integer, ney::Integer=-1,
    meshtype::Meshtype=QUADRANGLE; gmap=identity
)
    # Bounds
    i1, i2 = DomainSets.components(Ω)

    # Choose ny such that faces are nearly quadric if not specified
    if ney == -1
        w1 = IntervalSets.width(i1)
        w2 = IntervalSets.width(i2)
        ney = w2 / (w1 / nex) |> ceil |> Int
    end

    # Coordinates
    nnx = nex + 1
    nny = ney + 1
    nn = nnx * nny
    coordinates = zeros(2, nn)
    coordinates[1, :] = repeat(range(i1, length=nnx), outer=(ney + 1, 1))
    coordinates[2, :] = repeat(range(i2, length=nny), inner=(nex + 1, 1))

    # Apply geometric map
    if gmap != identity
        for i = 1:nn
            coordinates[:, i] = gmap(coordinates[:, i])
        end
    end

    # Geometry types
    g1 = g2 = GeometricObjectI
    if meshtype == QUADRANGLE && _isscaling(gmap)
        g2 = Box
    end

    # Mesh
    m = Mesh(coordinates, 2, g1=g1, g2=g2)

    # Connectivity
    cl = ConnectivityList()
    idx(i, j) = (j - 1) * nnx + i
    for j in 1:ney, i in 1:nex
        if meshtype == QUADRANGLE
            push!(cl, [idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)])
        elseif meshtype == TRIANGLE
            push!(cl, [idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)])
            push!(cl, [idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)])
        end
    end
    addlinks!(m.topology, 2, 0, cl)

    # Done
    return m
end

Meshes.Mesh(
    Ω::DomainSets.Rectangle, nex::Integer, ney::Integer=-1,
    meshtype::Meshtype=QUADRANGLE; gmap=identity
) = makemeshonrectangle(Ω, nex, ney, meshtype, gmap=gmap)

Meshes.Mesh(name::Symbol) = eval(:($name()))

function quadtri()
    coords = [0.0 1.0 2.0 0.1 0.9 1.9; 0.0 0.1 0.0 0.9 1.0 0.9]
    elts = [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
    return Mesh(coords, elts, 2)
end