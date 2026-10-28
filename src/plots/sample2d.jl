# Keyword arguments do not take part in dispatch, so dispatch on domain in positional methods
sample2d(f; domain, npoints::Integer, gmap=identity) = sample2d(f, domain, npoints, gmap)

sample2dlines(f; domain, npoints::Integer, nmeshlines=0, gmap=identity, zscale=1) =
    sample2dlines(f, domain, npoints, nmeshlines, gmap, zscale)

function sample2d(f, domain::DomainSets.Rectangle, npoints::Integer, gmap)
    n1, n2 = _getnintervals(domain, npoints)
    m = Mesh(domain, n1, n2, meshtype=TRIANGLE)
    x = tomatrix([[gmap(x)..., f(x)] for x in m.geometry.points])
    t = tomatrix([l for l in links(m.topology, 2, 0)], ROWS)
    return x, t
end

function sample2d(f, ::DomainSets.EuclideanUnitSimplex{2}, npoints::Integer, gmap)
    n = npoints

    # Lattice points (i/n, j/n) with i + j <= n, numbered row by row
    ps = [SA[i / n, j / n] for j = 0:n for i = 0:(n-j)]
    idx(i, j) = j * (n + 1) - j * (j - 1) ÷ 2 + i + 1

    # Two triangles per lattice cell, only one at the hypotenuse
    ts = Vector{Int}[]
    for j = 0:(n-1), i = 0:(n-1-j)
        push!(ts, [idx(i, j), idx(i + 1, j), idx(i, j + 1)])
        if i + j < n - 1
            push!(ts, [idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)])
        end
    end

    x = tomatrix([[gmap(p)..., f(p)] for p in ps])
    t = tomatrix(ts, ROWS)
    return x, t
end

function sample2dlines(f, domain::DomainSets.Rectangle, npoints::Integer, nmeshlines, gmap, zscale)
    points = []

    # Quick return
    isnothing(nmeshlines) && return points

    # Helpers
    ng1, ng2 = 2 .* _getnintervals(domain, npoints)
    nm1, nm2 = (nmeshlines isa Integer ? [nmeshlines, nmeshlines] : nmeshlines) .+ 2
    makepoints(ξ1, ξ2) = [
        [gmap([ξ1[i], ξ2[i]])..., zscale * f([ξ1[i], ξ2[i]])] for i = eachindex(ξ1)
    ]

    # Collect 
    for ξ in range(domain.a[2], domain.b[2], length=nm2)
        ξ1 = range(domain.a[1], domain.b[1], length=ng1)
        ξ2 = [ξ for _ in 1:ng1]
        push!(points, makepoints(ξ1, ξ2))
    end
    for ξ in range(domain.a[1], domain.b[1], length=nm1)
        ξ1 = [ξ for _ in 1:ng2]
        ξ2 = range(domain.a[2], domain.b[2], length=ng2)
        push!(points, makepoints(ξ1, ξ2))
    end

    return points
end

function sample2dlines(
    f, ::DomainSets.EuclideanUnitSimplex{2}, npoints::Integer, nmeshlines, gmap, zscale
)
    points = []

    # Quick return
    isnothing(nmeshlines) && return points

    # Helpers
    ng = 2 * npoints
    cs = range(0, 1, length=nmeshlines + 2)
    makeline(p1, p2) = [
        [gmap(p)..., zscale * f(p)] for p in (p1 + t * (p2 - p1) for t in range(0, 1, length=ng))
    ]

    # Lines parallel to the sides, the first or last value of c is a degenerate line
    for c in cs[1:end-1]
        push!(points, makeline(SA[0.0, c], SA[1 - c, c]))   # y = c
        push!(points, makeline(SA[c, 0.0], SA[c, 1 - c]))   # x = c
    end
    for c in cs[2:end]
        push!(points, makeline(SA[c, 0.0], SA[0.0, c]))     # x + y = c
    end

    return points
end
