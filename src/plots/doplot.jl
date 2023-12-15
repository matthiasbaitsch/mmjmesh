# -------------------------------------------------------------------------------------------------
# Main plot function
# -------------------------------------------------------------------------------------------------

"""
    xmplot(m::Mesh[, ps::PlotStyle=PlotStyle()])

Plot mesh `m` with style `ps`.
"""
function doplot(plot::MPlot, m::Mesh, ps::PlotStyle)

    show(s) = !isnothing(s) && s.visible

    if show(ps.lineplot)
        doplotlineplot(plot, m, ps.lineplot)
    end
    if show(ps.faces)
        doplotfaces(plot, m, ps.faces)
    end
    if show(ps.edges)
        doplotedges(plot, m, ps.edges, false)
    end
    if show(ps.featureedges)
        doplotedges(plot, m, ps.featureedges, true)
    end
    if show(ps.nodes)
        Makie.scatter!(
            plot,
            coordinates(m),
            color=ps.nodes.color,
            markersize=ps.nodes.size
        )
    end

    return plot
end

# -------------------------------------------------------------------------------------------------
# Details
# -------------------------------------------------------------------------------------------------

function doplotlineplot(plot::MPlot, m::Mesh, ps::LineplotStyle)

    # TODO complete rewrite using functions
    Nn = nentities(m.topology, 0)
    Ne = nentities(m.topology, 1)
    values = ps.values
    edges = links(m.topology, 1, 0)
    coords = coordinates(m)

    # Get values right
    dv = size(values)
    if length(dv) == 1
        if dv[1] == Nn
            values = mapreduce(permutedims, vcat, [values[l] for l in links(m.topology, 1, 0)])'
        elseif dv[1] == Ne
            values = [values values]'
        else
            @notimplemented
        end
        dv = size(values)
    end

    # Check preconditions
    @assert length(dv) == 2
    @assert dv[1] == 2
    @assert dv[2] == Ne

    # sz = size(boundingbox(m.geometry))
    minx = minimum(coordinates(m)[1, :])
    maxx = maximum(coordinates(m)[1, :])
    miny = minimum(coordinates(m)[2, :])
    maxy = maximum(coordinates(m)[2, :])
    lx = maxx - minx
    ly = maxy - miny
    sz = sqrt(lx * ly)
    if sz == 0
        sz = max(lx, ly)
    end

    vmx = maximum(abs, values)
    a = -ps.scale * sz / vmx

    c = Vector{Float64}()
    xf = Vector{Float64}()
    yf = Vector{Float64}()
    tf = Vector{Vector{Int}}()

    xe = Vector{Float64}()
    ye = Vector{Float64}()

    for e ∈ 1:Ne
        n1 = edges[e, 1]
        n2 = edges[e, 2]
        x1 = coords[:, n1]
        x2 = coords[:, n2]
        d = (x2 - x1) / norm(x2 - x1)
        n = [-d[2], d[1]]
        v1 = values[1, e]
        v2 = values[2, e]
        p1 = x1 + a * v1 * n
        p2 = x2 + a * v2 * n
        bi = length(xf) + 1

        # Face
        push!(c, v1, v1, v2, v2)
        push!(xf, x1[1], p1[1], x2[1], p2[1])
        push!(yf, x1[2], p1[2], x2[2], p2[2])

        # Edge
        push!(xe, x1[1], p1[1], p2[1], x2[1], NaN)
        push!(ye, x1[2], p1[2], p2[2], x2[2], NaN)

        # Sign change
        if v1 * v2 < 1e-10
            # v1 + (v2 - v1) * s = 0
            s = v1 / (v1 - v2)
            x3 = x1 + s * (x2 - x1)
            push!(xf, x3[1])
            push!(yf, x3[2])
            push!(c, 0.0)
            push!(tf, [bi, bi + 1, bi + 4])
            push!(tf, [bi + 4, bi + 2, bi + 3])
        else
            push!(tf, [bi, bi + 1, bi + 3])
            push!(tf, [bi, bi + 2, bi + 3])
        end
    end

    if !isnothing(ps.faces.color)
        c = ps.faces.color
    end

    if ps.faces.visible
        x = [xf yf]'
        tf = mapreduce(permutedims, vcat, tf)
        Makie.mesh!(plot, x, tf, color=c, colormap=ps.faces.colormap)
    end

    if ps.outlines.visible
        Makie.lines!(plot, xe, ye, linewidth=ps.outlines.linewidth, color=ps.outlines.color)
    end
end

function doplotfaces(plot::MPlot, m::Mesh, ps::MeshStyle)

    # Test if we have data
    haveData = ps.color isa Vector

    # Helpers
    coords = coordinates(m)
    Nn = nnodes(m)
    Nf = nfaces(m)

    # Collect triangles    
    if !haveData || length(ps.color) == Nn                           # No color or nodal color
        tf = Vector{Vector{Int}}()
        for l in links(m.topology, 2, 0)
            if length(l) == 3
                push!(tf, l)
            elseif length(l) == 4
                push!(tf, [l[1], l[2], l[3]])
                push!(tf, [l[1], l[3], l[4]])
            end
        end
        x = coords
        tf = mapreduce(permutedims, vcat, tf)
        c = ps.color        
    elseif length(ps.color) == Nf                                    # Element color
        cnt = 1
        tf = Vector{Vector{Int}}()
        c = Vector{Float64}()
        xx = Vector{Float64}()
        yy = Vector{Float64}()
        for (i, l) ∈ enumerate(links(m.topology, 2, 0))
            if length(l) == 3
                append!(c, ps.color[i] * ones(3))
                append!(xx, coords[1, l[[1, 2, 3]]])
                append!(yy, coords[2, l[[1, 2, 3]]])
                push!(tf, cnt:cnt+2)
                cnt += 3
            elseif length(l) == 4
                append!(c, ps.color[i] * ones(6))
                append!(xx, coords[1, l[[1, 2, 3, 1, 3, 4]]])
                append!(yy, coords[2, l[[1, 2, 3, 1, 3, 4]]])
                push!(tf, cnt:cnt+2)
                push!(tf, cnt+3:cnt+5)
                cnt += 6
            else
                @notimplemented
            end
        end
        x = [xx yy]'
        tf = mapreduce(permutedims, vcat, tf)
    else
        @notimplemented
    end

    # Plot
    Makie.mesh!(plot, x, tf, color=c, colormap=ps.colormap)
end

function doplotedges(plot::MPlot, m::Mesh, ps::LineStyle, featureedgesonly::Bool)

    # Coordinates and links from edges to faces
    coords = coordinates(m)
    l12 = featureedgesonly ? links(m.topology, 1, 2) : ConnectivityList()

    # Collect coordinates - currently only 2D
    xx = Vector{Float64}()
    yy = Vector{Float64}()
    for (i, l) in enumerate(links(m.topology, 1, 0))
        if !featureedgesonly || length(l12, i) == 1
            x1 = coords[:, l[1]]
            x2 = coords[:, l[2]]
            push!(xx, x1[1], x2[1], NaN)
            push!(yy, x1[2], x2[2], NaN)
        end
    end

    # Plot
    Makie.lines!(plot, xx, yy, linewidth=ps.linewidth, color=ps.color)
end
