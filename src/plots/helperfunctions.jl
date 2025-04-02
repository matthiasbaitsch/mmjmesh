"""
    _getnintervals(r, inintervals) -> (n1, n2)

Returns the number of intervals for each coordinate direction such that
the `Rectangle` parameter `r` is divided into approximate squares. The longer side
of `r` is divided into `nintervals` intervals.
"""
function _getnintervals(r::DomainSets.Rectangle, nintervals::Integer)
    @assert nintervals >= 1
    l1, l2 = IntervalSets.width.(DomainSets.components(r))
    h = max(l1, l2) / nintervals
    n1 = Int(ceil(l1 / h))
    n2 = Int(ceil(l2 / h))
    return n1, n2
end


"""
    _collectlines(points) -> (x1, x2 [, x3])

Collect point arrays

```{julia}
points = [[[p111, p112], [p121, p122], ...], [[p211, p212], [p221, p222], ...], ...]
```

into arrays 

```{julia}
x1 = [p111, p121, ..., NaN, p211, p221, ..., NaN, ...]
x2 = [p112, p122, ..., NaN, p212, p222, ..., NaN, ...]
```

which can be passed to `Makie`'s `lines` function.
"""
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


"""
    _mergemeshes(meshes) -> mesh

Merges the meshes

```{julia}
[(coords, triangles), (coords, triangles), ...]
```

into a single mesh

```{julia}
(coords, triangles)
```
"""
function _mergemeshes(meshes)
    xx = [Float32[], Float32[], Float32[]]
    tt = [Int[], Int[], Int[]]
    for c in meshes
        xf, tf = c
        pos = length(xx[1])
        for i = 1:3
            append!(xx[i], xf[i, :])
            append!(tt[i], pos .+ tf[:, i])
        end
    end
    return stack(xx, dims=1), stack(tt)
end


function _getcolor(x::Matrix, color, zscale)
    if typeof(color) == Int && 1 <= color <= 3
        return x[color, :] / zscale
    end
    return color
end


function _samplefaces(mesh, makef, nintervals, nmeshlines, zscale)
    cf = []
    cl1 = []
    cl2 = []

    for face = faces(mesh)
        f = makef(face)
        gmap = parametrization(geometry(face))
        d = domain(gmap)

        # Faces
        push!(
            cf,
            sample2d(f, domain=d, npoints=2 * nintervals, gmap=gmap)
        )

        # Edges
        append!(
            cl1,
            sample2dlines(f, domain=d, npoints=nintervals, nmeshlines=0, gmap=gmap, zscale=zscale)
        )

        # Meshlines
        append!(
            cl2,
            sample2dlines(f, domain=d, npoints=nintervals, nmeshlines=nmeshlines, gmap=gmap, zscale=zscale)
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

_coeffs(v1, v2) = [(v1 + v2) / 2, (v2 - v1) / 2]

function _tofunctions(mesh, values)
    values = _collectvalues(mesh, values)
    if length(size(values)) == 1
        return [Polynomial([values[i]], IHat) for i in eachindex(values)]
    else
        return [Polynomial(_coeffs(values[:, i]...), IHat) for i in axes(values, 2)]
    end
end

