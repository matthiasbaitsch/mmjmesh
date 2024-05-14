using DomainSets: Rectangle


function makeinitialmesh(r::Rectangle, npoints::Integer)
    d1 = r.b[1] - r.a[1]
    d2 = r.b[2] - r.a[2]
    fn = npoints / max(d1, d2)
    n1 = Int(ceil(fn * d1))
    n2 = Int(ceil(fn * d2))
    m = makemeshonrectangle(d1, d2, n1, n2, TRIANGLE)
    m.geometry.points.coordinates .+= r.a
    return m
end

function sample2d(f::FunctionRnToR{2}; domain=None, npoints::Integer=30)
    if isnothing(domain)
        domain = domain(f)
    end
    m = makeinitialmesh(domain, npoints)
    x = tomatrix([[x[1], x[2], f(x)] for x in m.geometry.points])
    t = tomatrix([l for l in links(m.topology, 2, 0)], ROWS)
    return x, t
end
