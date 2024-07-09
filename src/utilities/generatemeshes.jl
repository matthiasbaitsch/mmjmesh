
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
    makemeshonrectangle(w::Number, h::Number, nx::Int, ny::Int, meshtype::Meshtype)

Generate a quad mesh on the domain ``[0, w] \\times [0, h]`` with `nx` elements in the 
``x`` direction and `ny` elements in the ``y`` direction.
"""
function makemeshonrectangle(w::Number, h::Number, nx::Int, ny::Int, mt::Meshtype=QUADRANGLE)

    # Coordinates
    nn = (nx + 1) * (ny + 1)
    c = zeros(2, nn)
    c[1, :] = repeat(range(0, stop=w, length=nx + 1), outer=(ny + 1, 1))
    c[2, :] = repeat(range(0, stop=h, length=ny + 1), inner=(nx + 1, 1))

    # Geometry type
    g2 = mt == QUADRANGLE ? Box : GeometricObjectI

    # Mesh
    m = Mesh(c, 2, g2=g2)

    # Connectivity
    cl = ConnectivityList()
    idx(i, j) = (j - 1) * (nx + 1) + i
    for j in 1:ny, i in 1:nx
        if mt == QUADRANGLE
            push!(cl, [idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)])
        elseif mt == TRIANGLE
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

