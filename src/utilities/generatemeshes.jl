
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

    # Housekeeping
    populatepredfinedgroups!(m)

    # Done
    return m
end


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
    w::Real, h::Real, nex::Integer, ney::Integer=-1, meshtype::Meshtype=QUADRANGLE
) = makemeshonrectangle((0 .. w) × (0 .. h), nex, ney, meshtype)

function makemeshonrectangle(
    Ω::Rectangle, nex::Integer, ney::Integer=-1, meshtype::Meshtype=QUADRANGLE
)
    # Bounds
    i1, i2 = components(Ω)

    # Choose ny such that faces are nearly quadric if not specified
    if ney == -1
        w1 = width(i1)
        w2 = width(i2)
        ney = w2 / (w1 / nex) |> ceil |> Int
    end

    # Coordinates
    nnx = nex + 1
    nny = ney + 1
    nn = nnx * nny
    coordinates = zeros(2, nn)
    coordinates[1, :] = repeat(range(i1, length=nnx), outer=(ney + 1, 1))
    coordinates[2, :] = repeat(range(i2, length=nny), inner=(nex + 1, 1))

    # Geometry type for faces
    g2 = (meshtype == QUADRANGLE ? Box : GeometricObjectI)

    # Mesh
    m = Mesh(coordinates, 2, g2=g2)

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

    # Housekeeping
    populatepredfinedgroups!(m)

    # Done
    return m
end

