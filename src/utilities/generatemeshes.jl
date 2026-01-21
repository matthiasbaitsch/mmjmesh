_isscaling(gmap) = (gmap == identity)
_isscaling(m::AffineMapping) = m.A[1, 2] == 0 && m.A[2, 1] == 0


"""
    Meshtype

Types of meshes on structured grids.
"""
@enum Meshtype begin
    QUADRANGLE
    TRIANGLE
    CRISSCROSS
end


"""
    Mesh(name::Symbol)
    Mesh(I::Interval, n::Integer; gmap=identity)
    Mesh(l::Real, w::Real, nex::Integer, ney::Integer=-1; meshtype::Meshtype=QUADRANGLE, gmap=identity)
    Mesh(Ω::Rectangle, nex::Integer, ney::Integer=-1; meshtype::Meshtype=QUADRANGLE, gmap=identity)
"""
Meshes.Mesh(
    name::Symbol
) = eval(:($name()))

Meshes.Mesh(
    I::IntervalSets.AbstractInterval, n::Integer;
    gmap=t -> [t; 0.0]
) = makemeshoninterval(IntervalSets.leftendpoint(I), IntervalSets.rightendpoint(I), n, gmap)

Meshes.Mesh(
    l::Real, w::Real, nex::Integer, ney::Integer=-1;
    meshtype::Meshtype=QUADRANGLE, gmap=identity
) = makemeshonrectangle((0 .. l) × (0 .. w), nex, ney, meshtype, gmap)

Meshes.Mesh(
    Ω::DomainSets.Rectangle, nex::Integer, ney::Integer=-1;
    meshtype::Meshtype=QUADRANGLE, gmap=identity
) = makemeshonrectangle(Ω, nex, ney, meshtype, gmap)


# -------------------------------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------------------------------

function quadtri()
    coords = [0.0 1.0 2.0 0.1 0.9 1.9; 0.0 0.1 0.0 0.9 1.0 0.9]
    elts = [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
    return Mesh(coords, elts, 2)
end


function makemeshoninterval(a::Number, b::Number, n::Int, gmap)
    # Coordinates
    c = reduce(hcat, gmap.(range(a, stop=b, length=n + 1)))

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


function makemeshonrectangle(
    Ω::DomainSets.Rectangle, nex::Integer, ney::Integer, meshtype::Meshtype, gmap
)
    # Bounds
    ix, iy = DomainSets.components(Ω)
    wx = IntervalSets.width(ix)
    wy = IntervalSets.width(iy)

    # Choose ny such that faces are nearly quadric if not specified
    if ney == -1
        ney = wy / (wx / nex) |> ceil |> Int
    end

    # Coordinates
    nnx = nex + 1
    nny = ney + 1
    c = zeros(2, nnx * nny)
    c[1, :] = repeat(range(ix, length=nnx), outer=(ney + 1, 1))
    c[2, :] = repeat(range(iy, length=nny), inner=(nex + 1, 1))

    # Inner nodes for crisscross mesh
    if meshtype == CRISSCROSS
        c2 = zeros(2, nex * ney)
        dx = wx / nex
        dy = wy / ney
        ix2 = IntervalSets.leftendpoint(ix) + dx / 2 .. IntervalSets.rightendpoint(ix) - dx / 2
        iy2 = IntervalSets.leftendpoint(iy) + dy / 2 .. IntervalSets.rightendpoint(iy) - dy / 2
        c2[1, :] = repeat(range(ix2, length=nex), outer=(ney, 1))
        c2[2, :] = repeat(range(iy2, length=ney), inner=(nex, 1))
        c = hcat(c, c2)
    end

    # Apply geometric map
    if gmap != identity
        c = mapslices(gmap, c, dims=1)
    end

    # Geometry types
    g1 = g2 = GeometricObjectI
    if meshtype == QUADRANGLE && _isscaling(gmap)
        g2 = Box
    end

    # Mesh
    m = Mesh(c, 2, g1=g1, g2=g2)

    # Connectivity
    cl = ConnectivityList()
    idx(i, j) = (j - 1) * nnx + i
    idx2(i, j) = nnx * nny + (j - 1) * nex + i
    for j in 1:ney, i in 1:nex
        if meshtype == QUADRANGLE
            push!(cl, [idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)])
        elseif meshtype == TRIANGLE
            push!(cl, [idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)])
            push!(cl, [idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)])
        elseif meshtype == CRISSCROSS
            push!(cl, [idx(i, j), idx(i + 1, j), idx2(i, j)])
            push!(cl, [idx(i + 1, j), idx(i + 1, j + 1), idx2(i, j)])
            push!(cl, [idx(i + 1, j + 1), idx(i, j + 1), idx2(i, j)])
            push!(cl, [idx(i, j + 1), idx(i, j), idx2(i, j)])
        else
            error("Unknown meshtype: " + meshtype)
        end
    end
    addlinks!(m.topology, 2, 0, cl)

    # Groups for nodes on edges
    definegroup!(m, 0, :b1, 1:nnx)
    definegroup!(m, 0, :b2, nnx:nnx:(nnx*nny))
    definegroup!(m, 0, :b3, (nnx*(nny-1)+1):(nnx*nny))
    definegroup!(m, 0, :b4, 1:nnx:(nnx*nny))

    # Done
    return m
end

