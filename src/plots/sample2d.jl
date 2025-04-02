function sample2d(f; domain::DomainSets.Rectangle, npoints::Integer, gmap=identity)
    n1, n2 = _getnintervals(domain, npoints)
    m = Mesh(domain, n1, n2, TRIANGLE)
    x = tomatrix([[gmap(x)..., f(x)] for x in m.geometry.points])
    t = tomatrix([l for l in links(m.topology, 2, 0)], ROWS)
    return x, t
end

function sample2dlines(
    f; domain::DomainSets.Rectangle, npoints::Integer, nmeshlines=0, gmap=identity, zscale=1
)
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
