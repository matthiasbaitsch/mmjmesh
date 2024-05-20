using DomainSets: Rectangle
using MMJMesh.Utilities: makemeshonrectangle, TRIANGLE

function _getnpoints(r::Rectangle, npoints::Integer)
    d1 = r.b[1] - r.a[1]
    d2 = r.b[2] - r.a[2]
    fn = npoints / max(d1, d2)
    n1 = Int(ceil(fn * d1))
    n2 = Int(ceil(fn * d2))
    return n1, n2
end

function _append!(
    l1::Vector, l2::Vector, l3::Vector,
    x1::AbstractVector, x2::AbstractVector, x3::AbstractVector
)
    append!(l1, x1)
    append!(l2, x2)
    append!(l3, x3)
    push!(l1, NaN)
    push!(l2, NaN)
    push!(l3, NaN)
end


function makeinitialmesh(r::Rectangle, npoints::Integer)
    d1 = r.b[1] - r.a[1]
    d2 = r.b[2] - r.a[2]
    n1, n2 = _getnpoints(r, npoints)
    m = makemeshonrectangle(d1, d2, n1, n2, TRIANGLE)
    m.geometry.points.coordinates .+= r.a
    return m
end


function sample2d(f::FunctionRnToR{2}; domain, npoints::Integer)
    m = makeinitialmesh(domain, npoints)
    x = tomatrix([[x[1], x[2], f(x)] for x in m.geometry.points])
    t = tomatrix([l for l in links(m.topology, 2, 0)], ROWS)
    return x, t
end


function sample2dlines(d::Rectangle, f::FunctionRnToR{2}, npoints::Integer, mesh)
    l1s = Float64[]
    l2s = Float64[]
    l3s = Float64[]

    if isnothing(mesh)         
        return l1s, l2s, l3s
    end
    
    mesh = typeof(mesh) == Int ? [mesh, mesh] : mesh

    nm1, nm2 = mesh .+ 2
    ng1, ng2 = _getnpoints(d, npoints)

    ng1 *= 2
    ng2 *= 2

    for x1 in range(d.a[1], d.b[1], length=nm1)
        x1s = [x1 for _ in 1:ng2]
        x2s = range(d.a[2], d.b[2], length=ng2)
        x3s = f.([[x1s[i], x2s[i]] for i ∈ 1:ng2])
        _append!(l1s, l2s, l3s, x1s, x2s, x3s)
    end

    for x2 in range(d.a[2], d.b[2], length=nm2)
        x1s = range(d.a[1], d.b[1], length=ng1)
        x2s = [x2 for _ in 1:ng1]
        x3s = f.([[x1s[i], x2s[i]] for i ∈ 1:ng1])
        _append!(l1s, l2s, l3s, x1s, x2s, x3s)
    end

    return l1s, l2s, l3s
end