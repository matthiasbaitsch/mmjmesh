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
        featureedgecolor=MakieCore.automatic,
        featureedgelinewidth=MakieCore.automatic,

        # Faces
        facesvisible=MakieCore.automatic,
        facecolor=MakieCore.automatic,
        facecolormap=MakieCore.automatic,

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
        if colorbar && !scene.plot[:colorbygroups][]
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
    isautomatic(key) = (plot.attributes[key][] == MakieCore.automatic)
    setifautomatic(key, value) = isautomatic(key) && (plot.attributes[key] = value)

    # Flag if we have data
    color = length(plot) > 1 ? plot.scalars[] : plot.facecolor
    plot.havedata = color isa Vector

    # Configure
    plot.colorbygroups = false
    plot.lineplotvisible = false
    setifautomatic(:nodesvisible, nnodes(mesh) <= 50)

    if pdim(mesh) == 1
        plot.edgesvisible = true
        plot.edgecolor = MakieCore.theme(plot, :linecolor)
        plot.featureedgesvisible = false
        plot.facesvisible = false
        plot.lineplotvisible = length(plot) > 1
        setifautomatic(:edgelinewidth, 3)
    elseif pdim(mesh) == 2
        plot.colorbygroups = !plot.havedata[] && isautomatic(:facecolor) && hasgroups(mesh.groups, d=2)
        setifautomatic(:featureedgesvisible, true)
        setifautomatic(:edgesvisible, nfaces(mesh) <= 100)
        setifautomatic(:edgecolor, MakieCore.theme(plot, :linecolor))
        setifautomatic(:edgelinewidth, 0.75)
        setifautomatic(:featureedgelinewidth, 3 * plot.edgelinewidth[])
        setifautomatic(:facesvisible, true)
        if plot.colorbygroups[]
            setifautomatic(:facecolormap, :Pastel1_9)
        else
            setifautomatic(:facecolor, :seashell2)
            setifautomatic(:facecolormap, MakieCore.theme(plot, :colormap))
        end
    else
        @notimplemented
    end

    # Plot
    plot.lineplotvisible[] && plotlineplot(plot)
    plot.facesvisible[] && plotfaces(plot)
    plot.edgesvisible[] && plotedges(plot, false)
    plot.featureedgesvisible[] && plotedges(plot, true)
    plot.nodesvisible[] && MakieCore.scatter!(
        plot, coordinates(mesh), color=plot.nodecolor, markersize=plot.nodesize
    )

    # Return
    return plot
end


# -------------------------------------------------------------------------------------------------
# Plot parts
# -------------------------------------------------------------------------------------------------

plotlineplot(plot::MPlot) =
    plotlineplot(plot, plot.mesh[], plot.scalars[])

# Version for numbers
plotlineplot(plot::MPlot, mesh::Mesh, values::AbstractVecOrMat{<:Real}) =
    plotlineplot(plot, mesh, _tofunctions(mesh, values))

# Version for functions
function plotlineplot(plot::MPlot, mesh::Mesh, functions::AbstractVector{<:FunctionRToR{IHat}})

    # Check input
    @assert nedges(mesh) == length(functions)

    # Scaling factor
    fmax(f) = maximum(abs.(f.(-1:1:1)))
    vmx = maximum(fmax.(functions))
    sz = diagonal(boundingbox(mesh.geometry))
    a = -plot.lineplotscale[] * sz / vmx

    # Storage
    xe = Float32[]
    ye = Float32[]
    cf = Float32[]
    xf = Float32[]
    yf = Float32[]
    triangles = Int[]

    # Loop over edges
    for i ∈ eachindex(functions)

        # Mappings for edge, function and line plot
        f = functions[i]
        ce = parametrization(geometry(edge(mesh, i)))
        cl = ce + a * f * UnitNormal(ce)

        # Sample, refactor for curved edges
        params, values = sampleadaptive(f, -1.0, 1.0, ir=true, rp=true, yscale=abs(a))
        pe = tomatrix(ce.(params))
        pl = tomatrix(cl.(params))

        # Append
        _appendedges!(xe, ye, pe, pl)
        _appendfaces!(xf, yf, cf, triangles, pe, pl, values[2, :])
    end

    # Color for faces if specified
    if plot.lineplotfacescolor[] != MakieCore.automatic
        cf = plot.lineplotfacescolor[]
    end

    # Plot faces
    if plot.lineplotfacesvisible[]
        MakieCore.mesh!(plot, [xf yf]', reshape(triangles, 3, :)',
            color=cf, colormap=plot.lineplotfacescolormap)
    end

    # Plot outlines
    if plot.lineplotoutlinesvisible[]
        MakieCore.lines!(plot, xe, ye, linewidth=plot.lineplotoutlineslinewidth,
            color=plot.lineplotoutlinescolor)
    end
end


function plotfaces(plot::MPlot)
    mesh = plot.mesh[]
    color = length(plot) > 1 ? plot.scalars[] : plot.facecolor

    # Test if we have data and overall color
    haveData = color isa Vector
    haveColor = plot.facecolor[] != MakieCore.automatic

    # Helpers
    coords = coordinates(mesh)
    Nn = nnodes(mesh)
    Nf = nfaces(mesh)

    # Color by groups
    if !haveData && !haveColor && hasgroups(mesh.groups, d=2)
        color = groupids(mesh, d=2, predefined=false)
        haveData = true
    end

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
    MakieCore.mesh!(plot, x, tf, color=c, colormap=plot.facecolormap)
end


function plotedges(plot::MPlot, featureedges::Bool)
    mesh = plot.mesh[]

    # Collect indices of edges to plot
    if featureedges
        indices = mesh.groups[:boundaryedges]
        for n in groupnames(mesh.groups, d=1, predefined=false)
            indices = indices ∪ mesh.groups[n]
        end
    else
        indices = 1:nedges(mesh)
    end

    # Collect coordinates
    # TODO 3D: Generalize
    xx = Real[]
    yy = Real[]
    for i ∈ indices
        e = edge(mesh, i)
        x1 = coordinates(e, 1)
        x2 = coordinates(e, 2)
        push!(xx, x1[1], x2[1], NaN)
        push!(yy, x1[2], x2[2], NaN)
    end

    # Plot attributes
    if featureedges
        lc = plot.featureedgecolor[]
        lw = plot.featureedgelinewidth[]
        if lc == MakieCore.automatic
            if ngroups(mesh.groups, d=1) > 0
                ids = groupids(mesh, d=1, predefined=true)
                lc = reshape(repeat(ids[indices], 1, 3)', :)
            else
                lc = :black
            end
        end
    else
        lw = plot.edgelinewidth[]
        lc = plot.edgecolor[]
    end

    # Plot
    MakieCore.lines!(plot, xx, yy, linewidth=lw, color=lc, colormap=:tab10)
end


# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

function _appendedges!(xe, ye, edgepoints, lineplotpoints)
    push!(xe, edgepoints[1, 1])
    append!(xe, lineplotpoints[1, :])
    append!(xe, [edgepoints[1, end], NaN])
    push!(ye, edgepoints[2, 1])
    append!(ye, lineplotpoints[2, :])
    append!(ye, [edgepoints[2, end], NaN])
end

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

