# -------------------------------------------------------------------------------------------------
# Recipe
# -------------------------------------------------------------------------------------------------

MakieCore.@recipe(MPlot, mesh, scalars) do scene

    # Collect attributes from theme
    attr = MakieCore.Attributes(

        # Nodes
        nodesvisible=MakieCore.theme(scene, :nodesvisible),
        nodecolor=MakieCore.theme(scene, :nodecolor),
        nodesize=MakieCore.theme(scene, :markersize),

        # Edges
        edgesvisible=MakieCore.theme(scene, :edgesvisible),
        edgecolor=MakieCore.theme(scene, :edgecolor),
        edgelinewidth=MakieCore.theme(scene, :edgelinewidth),

        # Featureedges
        featureedgesvisible=MakieCore.theme(scene, :featureedgesvisible),
        featureedgecolor=MakieCore.theme(scene, :featureedgecolor),
        featureedgelinewidth=MakieCore.theme(scene, :featureedgelinewidth),

        # Faces
        facesvisible=MakieCore.theme(scene, :facesvisible),
        facecolor=MakieCore.theme(scene, :facecolor),
        facecolormap=MakieCore.theme(scene, :facecolormap),

        # Lineplot
        lineplotvisible=MakieCore.theme(scene, :lineplotvisible),
        lineplotvalues=MakieCore.theme(scene, :lineplotvalues),
        lineplotscale=MakieCore.theme(scene, :lineplotscale),
        lineplotoutlinesvisible=MakieCore.theme(scene, :lineplotoutlinesvisible),
        lineplotoutlinescolor=MakieCore.theme(scene, :lineplotoutlinescolor),
        lineplotoutlineslinewidth=MakieCore.theme(scene, :lineplotoutlineslinewidth),
        lineplotfacesvisible=MakieCore.theme(scene, :lineplotfacesvisible),
        lineplotfacescolor=MakieCore.theme(scene, :lineplotfacescolor),
        lineplotfacescolormap=MakieCore.theme(scene, :lineplotfacescolormap),

        # Faceplot
        faceplotzscale=MakieCore.theme(scene, :faceplotzscale),
        faceplotmesh=MakieCore.theme(scene, :faceplotmesh),
        faceplotnpoints=MakieCore.theme(scene, :faceplotnpoints),
        faceplotmeshcolor=MakieCore.theme(scene, :faceplotmeshcolor),
        faceplotmeshlinewidth=MakieCore.theme(scene, :faceplotmeshlinewidth),
    )

    MakieCore.generic_plot_attributes!(attr)
    MakieCore.colormap_attributes!(attr, MakieCore.theme(scene, :colormap))

    return attr
end


# -------------------------------------------------------------------------------------------------
# Main plot function
# -------------------------------------------------------------------------------------------------

function MakieCore.plot!(plot::MPlot)
    mesh = plot.mesh[]

    # Helpers
    isundefined(key) = key ∉ keys(plot.attributes) || isnothing(plot.attributes[key][])
    setifundefined(key, value) = isundefined(key) && (plot.attributes[key] = value)

    # Configure
    color = length(plot) > 1 ? plot.scalars[] : plot.facecolor
    colorbygroups = false
    lineplotvisible = false

    # Settings if not defined
    setifundefined(:nodesvisible, nnodes(mesh) <= 50)
    setifundefined(:nodecolor, :tomato)
    setifundefined(:edgecolor, MakieCore.theme(plot, :linecolor))
    setifundefined(:lineplotscale, 0.1)
    setifundefined(:lineplotoutlinesvisible, true)
    setifundefined(:lineplotfacesvisible, true)
    setifundefined(:lineplotoutlinescolor, MakieCore.theme(plot, :linecolor))
    setifundefined(:lineplotoutlineslinewidth, MakieCore.theme(plot, :linewidth))
    setifundefined(:lineplotfacescolormap, MakieCore.theme(plot, :colormap))
    setifundefined(:faceplotzscale, 0.0)
    setifundefined(:faceplotnpoints, 30)
    setifundefined(:faceplotmeshlinewidth, 1.25)
    setifundefined(:faceplotmeshcolor, MakieCore.theme(plot, :linecolor))

    # Access node coordinates or warp
    if isundefined(:nodewarp)
        plot.nodecoordinates = node -> coordinates(node)
    else
        nodewarp = plot.nodewarp[]
        if nodewarp isa Function
            plot.nodecoordinates = plot.nodewarp[]
        elseif nodewarp isa AbstractArray
            plot.nodecoordinates = (node -> [coordinates(node)..., nodewarp[index(node)]])
        end
    end

    # Settings for mesh of 1D elements
    if pdim(mesh) == 1
        plot.edgesvisible = true
        plot.featureedgesvisible = false
        plot.facesvisible = false
        lineplotvisible = length(plot) > 1
        setifundefined(:edgelinewidth, 3)
    end

    # Settings for mesh of 2D elements
    if pdim(mesh) == 2
        havedata = color isa Vector
        colorbygroups = !havedata && isundefined(:facecolor) && hasgroups(mesh.groups, d=2)
        setifundefined(:featureedgesvisible, true)
        setifundefined(:edgesvisible, nfaces(mesh) <= 100)
        setifundefined(:edgecolor, MakieCore.theme(plot, :linecolor))
        setifundefined(:edgelinewidth, 0.75)
        setifundefined(:featureedgelinewidth, 3 * plot.edgelinewidth[])
        setifundefined(:facesvisible, true)
        if colorbygroups
            setifundefined(:facecolormap, :Pastel1_9)
        else
            if color isa Function || color isa Symbol
                setifundefined(:facecolor, 3)
            else
                setifundefined(:facecolor, :seashell2)
            end
            setifundefined(:facecolormap, MakieCore.theme(plot, :colormap))
        end
    end

    # Higher dimensions not implemented yet
    if pdim(mesh) >= 3
        @notimplemented
    end

    # Used by mconf
    plot.colorbygroups = colorbygroups

    # Plot functions on faces goes extra at the moment - TODO this is a hack, refactor
    if pdim(mesh) == 2 && (color isa Function || color isa Symbol)
        plotfacefunctions(plot)
    else # Plot
        lineplotvisible && plotlineplot(plot)
        plot.facesvisible[] && plotfaces(plot)
        plot.edgesvisible[] && plotedges(plot, false)
        plot.featureedgesvisible[] && plotedges(plot, true)
        plot.nodesvisible[] && plotnodes(plot)
    end

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
        params, values = sample1d(f, -1.0, 1.0, ir=true, rp=true, yscale=abs(a))
        pe = tomatrix(ce.(params))
        pl = tomatrix(cl.(params))

        # Append
        _appendedges!(xe, ye, pe, pl)
        _appendfaces!(xf, yf, cf, triangles, pe, pl, values[2, :])
    end

    # Color for faces if specified
    if !isnothing(plot.lineplotfacescolor[])
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

    # Mesh and color
    mesh = plot.mesh[]
    nodecoordinates = plot.nodecoordinates[]
    color = length(plot) > 1 ? plot.scalars[] : plot.facecolor

    # Test if we have data and overall color
    haveData = color isa Vector
    haveColor = !isnothing(plot.facecolor[])

    # Helpers
    coords = nodecoordinates.(nodes(mesh))
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
    elseif length(color) == Nf                                    # Element color
        cnt = 1
        tf = []
        coordstmp = []
        colortmp = Float32[]
        for (i, l) ∈ enumerate(links(mesh.topology, 2, 0))
            if length(l) == 3
                append!(colortmp, color[i] * ones(3))
                append!(coordstmp, coords[l[[1, 2, 3]]])
                push!(tf, cnt:cnt+2)
                cnt += 3
            elseif length(l) == 4
                append!(colortmp, color[i] * ones(6))
                append!(coordstmp, coords[l[[1, 2, 3, 1, 3, 4]]])
                push!(tf, cnt:cnt+2)
                push!(tf, cnt+3:cnt+5)
                cnt += 6
            else
                @notimplemented
            end
        end
        color = colortmp
        coords = coordstmp
    else
        @notimplemented
    end

    # Plot
    MakieCore.mesh!(
        plot,
        tomatrix(coords), tomatrix(tf, ROWS),
        color=color, colormap=plot.facecolormap, colorrange=plot.colorrange
    )
end


function plotedges(plot::MPlot, featureedges::Bool)
    mesh = plot.mesh[]
    nodecoordinates = plot.nodecoordinates[]

    # Collect indices of edges to plot
    if featureedges
        indices = mesh.groups[:boundaryedges]
        for n in groupnames(mesh.groups, d=1, predefined=false)
            indices = indices ∪ mesh.groups[n]
        end
    else
        indices = 1:nedges(mesh)
    end

    # Coordinates
    points = [nodecoordinates.(nodes(edge(mesh, i))) for i = indices]

    # Plot attributes
    if featureedges
        lc = plot.featureedgecolor[]
        lw = plot.featureedgelinewidth[]
        if isnothing(lc)
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

    MakieCore.lines!(plot, _collectlines(points)..., linewidth=lw, color=lc, colormap=:tab10)
end

function plotnodes(plot::MPlot)
    mesh = plot.args[1][]
    nodecoordinates = plot.nodecoordinates[]
    x = nodecoordinates.(nodes(mesh)) |> tomatrix

    MakieCore.scatter!(
        plot, x, color=plot.nodecolor, markersize=plot.nodesize
    )
end

function plotfacefunctions(plot::MPlot)
    attributes = plot.attributes

    # Edge attributes
    edgesvisible = attributes.edgesvisible[]
    edgecolor = attributes.edgecolor[]
    edgelinewidth = attributes.edgelinewidth[]

    # Face attributes
    facecolor = attributes.facecolor[]

    # Faceplot
    fzscale = attributes.faceplotzscale[]
    fnpoints = attributes.faceplotnpoints[]
    fmesh = attributes.faceplotmesh[]
    fmeshcolor = attributes.faceplotmeshcolor[]
    fmeshlinewidth = attributes.faceplotmeshlinewidth[]

    # General
    colorrange = attributes.colorrange[]
    colormap = attributes.colormap[]

    # Parameters
    mesh = plot.args[1][]
    colors = plot.args[2][]

    # Function to plot
    if colors isa Symbol
        facefunction(face) = face.data[:post](face, colors)
    else
        facefunction = colors
    end

    # Collect
    cf, clEdges, clMesh = _samplefaces(mesh, facefunction, fnpoints, fzscale, fmesh)
    xf, tf = _collectfaces(cf)
    color = _getcolor(xf, facecolor, 1)
    xf[3, :] *= fzscale

    # Plot Mesh
    mesh!(
        plot, xf, tf, color=color,
        colormap=colormap, colorrange=colorrange
    )

    # Plot mesh lines on elements
    lines!(plot, _collectlines(clMesh)..., color=fmeshcolor, linewidth=fmeshlinewidth)

    # Plot element edges
    if edgesvisible
        lines!(plot, _collectlines(clEdges)..., color=edgecolor, linewidth=edgelinewidth)
    end
end


# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

# TODO move to helper functions file

function _samplefaces(m, mf, npoints, zscale, nmeshlines)
    cf = []
    cl1 = []
    cl2 = []
    for face = faces(m)
        f = mf(face)
        gmap = _makegmap(face)
        push!(cf, sample2d(f, domain=QHat, npoints=2 * npoints, gmap=gmap))
        append!(cl1, sample2dlines(f, domain=QHat, npoints=npoints, nmeshlines=0, gmap=gmap, zscale=zscale))
        append!(cl2, sample2dlines(f, domain=QHat, npoints=npoints, nmeshlines=nmeshlines, gmap=gmap, zscale=zscale))
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

function _fsize(face)
    x = coordinates(face)
    p = x[:, 1]
    l1 = x[1, 2] - x[1, 1]
    l2 = x[2, 3] - x[2, 2]
    return p, l1, l2
end

function _makegmap(face) # TODO use geometry map for face
    nn(x) = (1 + x) / 2
    p, a, b = _fsize(face)
    return x -> p + [nn(x[1]) * a, nn(x[2]) * b]
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


