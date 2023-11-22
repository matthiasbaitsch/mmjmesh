"""
    makemeshonrectangle(w::Number, h::Number, nx::Int, ny::Int)

Generate a quad mesh on the domain ``[0, w] \\times [0, h]`` with `nx` elements in the 
``x`` direction and `ny` elements in the ``y`` direction.
"""
function makemeshonrectangle(w::Number, h::Number, nx::Int, ny::Int)

    # Coordinates
    nn = (nx + 1) * (ny + 1)
    c = zeros(2, nn)
    c[1, :] = repeat(range(0, stop=w, length=nx + 1), outer=(ny + 1, 1))
    c[2, :] = repeat(range(0, stop=h, length=ny + 1), inner=(nx + 1, 1))

    # Mesh
    m = Mesh(c)

    # Connectivity
    cl = ConnectivityList()
    idx(i, j) = (j - 1) * (nx + 1) + i
    for j in 1:ny, i in 1:nx
        push!(cl, [idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)])
    end
    addlinks!(m.topology, 2, 0, cl)

    # Done
    return m
end

