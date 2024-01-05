# -------------------------------------------------------------------------------------------------
# Recipe
# -------------------------------------------------------------------------------------------------

MakieCore.@recipe(MPlot, mesh, scalars) do scene

    attr = MakieCore.Attributes(

        # Nodes
        nodesvisible=MakieCore.automatic,
        nodecolor=:tomato,
        nodesize=MakieCore.theme(scene, :markersize),

        # Edges
        edgesvisible=MakieCore.automatic,
        edgecolor=MakieCore.automatic,
        edgelinewidth=MakieCore.automatic,

        # Featureedges
        featureedgesvisible=MakieCore.automatic,
        featureedgecolor=MakieCore.theme(scene, :linecolor),
        featureedgelinewidth=MakieCore.automatic,

        # Faces
        facesvisible=MakieCore.automatic,
        facecolor=:seashell2,
        facecolormap=MakieCore.theme(scene, :colormap),

        # Lineplot
        lineplotvisible=MakieCore.automatic,
        lineplotvalues=MakieCore.automatic,
        lineplotscale=0.1,
        lineplotoutlinesvisible=true,
        lineplotoutlinescolor=MakieCore.theme(scene, :linecolor),
        lineplotoutlineslinewidth=MakieCore.theme(scene, :linewidth),
        lineplotfacesvisible=true,
        lineplotfacescolor=MakieCore.automatic,
        lineplotfacescolormap=MakieCore.theme(scene, :colormap),
    )

    MakieCore.generic_plot_attributes!(attr)
    MakieCore.colormap_attributes!(attr, MakieCore.theme(scene, :colormap))

    return attr
end

# -------------------------------------------------------------------------------------------------
# Configure
# -------------------------------------------------------------------------------------------------

function mconf(; colorbar=true, dataaspect=true, blank=true, title="")
    return scene -> begin
        ax = scene.axis
        if dataaspect
            ax.aspect = Makie.DataAspect()
        end
        if blank
            Makie.hidedecorations!(ax)
            Makie.hidespines!(ax)
        end
        if colorbar
            # Hack using the fact that the mesh is always plotted first
            p = scene.plot.plots[1]
            if !isnothing(Makie.extract_colormap_recursive(p))
                Makie.Colorbar(scene.figure[1, 2], p)
            end
        end
        if title != ""
            ax.title = title
        end
        return scene
    end
end

# -------------------------------------------------------------------------------------------------
# Main plot function
# -------------------------------------------------------------------------------------------------

function MakieCore.plot!(plot::MPlot)
    mesh = plot.mesh[]

    # Helper
    setifautomatic(key, value) = (plot.attributes[key][] == MakieCore.automatic) && (plot.attributes[key] = value)

    # Configure
    setifautomatic(:nodesvisible, nnodes(mesh) <= 50)
    if pdim(mesh) == 1                              # 1D mesh
        plot.edgesvisible = true
        plot.edgecolor = MakieCore.theme(plot, :linecolor)
        plot.featureedgesvisible = false
        plot.facesvisible = false
        if length(plot) > 1
            plot.lineplotvisible = true
        else
            plot.lineplotvisible = false
        end
        setifautomatic(:edgelinewidth, 3)
    elseif pdim(mesh) == 2                          # 2D mesh
        plot.lineplotvisible = false
        setifautomatic(:featureedgesvisible, true)
        setifautomatic(:edgesvisible, nfaces(mesh) <= 100)
        setifautomatic(:edgecolor, :gray50)
        setifautomatic(:edgelinewidth, 0.75)
        setifautomatic(:featureedgelinewidth, 3 * plot.edgelinewidth[])
        setifautomatic(:facesvisible, true)
    else
        @notimplemented
    end

    # Plot
    if plot.lineplotvisible[]
        plotlineplot(plot)
    end
    if plot.facesvisible[]
        plotfaces(plot)
    end
    if plot.edgesvisible[]
        plotedges(plot, false)
    end
    if plot.featureedgesvisible[]
        plotedges(plot, true)
    end
    if plot.nodesvisible[]
        Makie.scatter!(
            plot,
            coordinates(mesh),
            color=plot.nodecolor,
            markersize=plot.nodesize
        )
    end

    # Return
    return plot
end

# -------------------------------------------------------------------------------------------------
# Details
# -------------------------------------------------------------------------------------------------

function plotlineplot(plot::MPlot)
    mesh = plot.mesh[]

    # TODO complete rewrite using functions
    Nn = nentities(mesh.topology, 0)
    Ne = nentities(mesh.topology, 1)
    values = plot.scalars[]
    edges = links(mesh.topology, 1, 0)
    coords = coordinates(mesh)

    # Get values right
    dv = size(values)
    if length(dv) == 1
        if dv[1] == Nn
            values = mapreduce(permutedims, vcat, [values[l] for l in links(mesh.topology, 1, 0)])'
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
    minx = minimum(coordinates(mesh)[1, :])
    maxx = maximum(coordinates(mesh)[1, :])
    miny = minimum(coordinates(mesh)[2, :])
    maxy = maximum(coordinates(mesh)[2, :])
    lx = maxx - minx
    ly = maxy - miny
    sz = sqrt(lx * ly)
    if sz == 0
        sz = max(lx, ly)
    end

    vmx = maximum(abs, values)
    a = -plot.lineplotscale[] * sz / vmx

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

    if plot.lineplotfacescolor[] != MakieCore.automatic
        c = plot.lineplotfacescolor[]
    end

    if plot.lineplotfacesvisible[]
        x = [xf yf]'
        tf = mapreduce(permutedims, vcat, tf)
        Makie.mesh!(plot, x, tf, color=c, colormap=plot.lineplotfacescolormap)
    end

    if plot.lineplotoutlinesvisible[]
        Makie.lines!(plot, xe, ye, linewidth=plot.lineplotoutlineslinewidth, color=plot.lineplotoutlinescolor)
    end
end

function plotfaces(plot::MPlot)
    mesh = plot.mesh[]
    color = length(plot) > 1 ? plot.scalars[] : plot.facecolor

    # Test if we have data
    haveData = color isa Vector

    # Helpers
    coords = coordinates(mesh)
    Nn = nnodes(mesh)
    Nf = nfaces(mesh)

    # Collect triangles    
    if !haveData || length(color) == Nn                           # No color or nodal color
        tf = Vector{Vector{Int}}()
        for l in links(mesh.topology, 2, 0)
            if length(l) == 3
                push!(tf, l)
            elseif length(l) == 4
                push!(tf, [l[1], l[2], l[3]])
                push!(tf, [l[1], l[3], l[4]])
            end
        end
        x = coords
        tf = mapreduce(permutedims, vcat, tf)
        c = color
    elseif length(color) == Nf                                    # Element color
        cnt = 1
        tf = Vector{Vector{Int}}()
        c = Vector{Float64}()
        xx = Vector{Float64}()
        yy = Vector{Float64}()
        for (i, l) ∈ enumerate(links(mesh.topology, 2, 0))
            if length(l) == 3
                append!(c, color[i] * ones(3))
                append!(xx, coords[1, l[[1, 2, 3]]])
                append!(yy, coords[2, l[[1, 2, 3]]])
                push!(tf, cnt:cnt+2)
                cnt += 3
            elseif length(l) == 4
                append!(c, color[i] * ones(6))
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
    Makie.mesh!(plot, x, tf, color=c, colormap=plot.facecolormap)
end

function plotedges(plot::MPlot, featureedges::Bool)
    mesh = plot.mesh[]

    # Coordinates and links from edges to faces
    coords = coordinates(mesh)
    l12 = featureedges ? links(mesh.topology, 1, 2) : ConnectivityList()

    # Collect coordinates - currently only 2D
    xx = Vector{Float64}()
    yy = Vector{Float64}()
    for (i, l) in enumerate(links(mesh.topology, 1, 0))
        if !featureedges || length(l12, i) == 1
            x1 = coords[:, l[1]]
            x2 = coords[:, l[2]]
            push!(xx, x1[1], x2[1], NaN)
            push!(yy, x1[2], x2[2], NaN)
        end
    end

    # Plot attributes
    if featureedges
        lw = plot.featureedgelinewidth[]
        lc = plot.featureedgecolor[]
    else
        lw = plot.edgelinewidth[]
        lc = plot.edgecolor[]
    end

    # Plot
    Makie.lines!(plot, xx, yy, linewidth=lw, color=lc)
end