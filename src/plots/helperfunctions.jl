function _getnpoints(r::Rectangle, npoints::Integer)
    d1 = r.b[1] - r.a[1]
    d2 = r.b[2] - r.a[2]
    fn = npoints / max(d1, d2)
    n1 = Int(ceil(fn * d1))
    n2 = Int(ceil(fn * d2))
    return n1, n2
end


function _samplefaces(m, mf, npoints, zscale, nmeshlines)
    cf = []
    cl1 = []
    cl2 = []

    for face = faces(m)
        f = mf(face)
        gmap = parametrization(geometry(face))
        d = domain(gmap)

        # Faces
        push!(
            cf,
            sample2d(f, domain=d, npoints=2 * npoints, gmap=gmap)
        )

        # Edges
        append!(
            cl1,
            sample2dlines(f, domain=d, npoints=npoints, nmeshlines=0, gmap=gmap, zscale=zscale)
        )

        # Meshlines
        append!(
            cl2,
            sample2dlines(f, domain=d, npoints=npoints, nmeshlines=nmeshlines, gmap=gmap, zscale=zscale)
        )
    end

    return cf, cl1, cl2
end

# TODO refactor to get rid of this
function _appendedges!(xe, ye, edgepoints, lineplotpoints)
    push!(xe, edgepoints[1, 1])
    append!(xe, lineplotpoints[1, :])
    append!(xe, [edgepoints[1, end], NaN])
    push!(ye, edgepoints[2, 1])
    append!(ye, lineplotpoints[2, :])
    append!(ye, [edgepoints[2, end], NaN])
end

# TODO refactor to get rid of this
function _appendfaces!(xf, yf, cf, triangles, edgepoints, lineplotpoints, values)
    bi = length(xf) + 1
    np = length(values)
    append!(xf, edgepoints[1, :])
    append!(xf, lineplotpoints[1, :])
    append!(yf, edgepoints[2, :])
    append!(yf, lineplotpoints[2, :])
    append!(cf, values)
    append!(cf, values)
    for j ∈ 1:np-1
        append!(triangles, [bi + j - 1, bi + j, bi + np + j - 1])
        append!(triangles, [bi + j, bi + np + j - 1, bi + np + j])
    end
end

function _collectvalues(mesh::Mesh, values)
    dv = size(values)
    Nn = nentities(mesh.topology, 0)
    Ne = nentities(mesh.topology, 1)

    # Check input
    @assert (length(dv) == 1 && (dv[1] == Nn || dv[1] == Ne)) ||
            (length(dv) == 2 && dv[1] == 2 && dv[2] == Ne)

    # Values from nodes to edges
    if length(dv) == 1 && dv[1] == nnodes(mesh)
        values = tomatrix([values[l] for l in links(mesh.topology, 1, 0)])
    end

    return values
end

function _collectlines(points)
    isempty(points) && return Float32[], Float32[]
    n = length(points[1][1])
    x1 = Float32[]
    x2 = Float32[]
    x3 = Float32[]
    for lp = points
        for p = lp
            push!(x1, p[1])
            push!(x2, p[2])
            n == 3 && push!(x3, p[3])
        end
        push!(x1, NaN)
        push!(x2, NaN)
        n == 3 && push!(x3, NaN)
    end
    n == 2 && return x1, x2
    n == 3 && return x1, x2, x3
end

function _collectfaces(cf)
    xx = [Float32[], Float32[], Float32[]]
    tt = [Int[], Int[], Int[]]
    for c in cf
        xf, tf = c
        pos = length(xx[1])
        for i = 1:3
            append!(xx[i], xf[i, :])
            append!(tt[i], pos .+ tf[:, i])
        end
    end
    return stack(xx, dims=1), stack(tt)
end

_coeffs(v1, v2) = [(v1 + v2) / 2, (v2 - v1) / 2]

function _tofunctions(mesh, values)
    values = _collectvalues(mesh, values)
    if length(size(values)) == 1
        return [Polynomial([values[i]], IHat) for i in eachindex(values)]
    else
        return [Polynomial(_coeffs(values[:, i]...), IHat) for i in axes(values, 2)]
    end
end

