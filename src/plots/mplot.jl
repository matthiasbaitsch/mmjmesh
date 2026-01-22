const DEFAULT_EDGE_LINEWIDTH_2D = 0.75
const DEFAULT_EDGE_COLORMAP = :tab10
const DEFAULT_FEATUREEDGE_LINEWIDTH_2D = 2.25
const DEFAULT_FACE_COLORCOMPONENT = 3
const DEFAULT_FACE_COLOR = :seashell2
const DEFAULT_FACE_COLORMAP_GROUPS = :Pastel1_9


# -------------------------------------------------------------------------------------------------
# Recipe
# -------------------------------------------------------------------------------------------------

Makie.@recipe MPlot (mesh, scalars_raw, vectors) begin

    # Nodes
    "Nodes visible."
    nodesvisible = automatic
    "Node color."
    nodecolor = :tomato
    "Node size."
    nodesize = @inherit markersize
    "Function or array to move nodes TODO."
    nodewarp = automatic

    # Edges
    "Edges visible."
    edgesvisible = automatic
    "Edge color. Use `:groups` to color by groups."
    edgecolor = automatic
    "Edge linewidth."
    edgelinewidth = automatic

    # Featureedges
    "Feature edges visible"
    featureedgesvisible = true
    "Feature edge color. Use `:groups` to color by groups."
    featureedgecolor = automatic
    "Feature edge linewidth."
    featureedgelinewidth = automatic

    # Lineplot
    lineplotvisible = automatic
    lineplotscale = 0.1
    lineplotoutlinesvisible = true
    lineplotoutlinescolor = :black
    lineplotoutlineslinewidth = @inherit linewidth
    lineplotfacesvisible = true
    lineplotfacescolor = automatic
    lineplotfacescolormap = @inherit colormap

    # Faces
    facesvisible = true
    facecolor = automatic
    facecolormap = automatic
    faceplotfunction = false

    # Faceplot
    faceplotzscale = 0.0
    faceplotmesh = 15
    faceplotnpoints = 30
    faceplotmeshcolor = @inherit linecolor
    faceplotmeshlinewidth = 1.25

    # Mist
    uscale = automatic
    smooth = false

    # Groups 
    "Use groups for coloring."
    groupsforcolors = true

    # Defaults
    Makie.mixin_colormap_attributes()...
    Makie.mixin_generic_plot_attributes()...
end

# -------------------------------------------------------------------------------------------------
# Main plot function
# -------------------------------------------------------------------------------------------------

Makie.convert_arguments(::Type{MPlot}, m::Mesh) = (m, nothing, nothing)

Makie.convert_arguments(::Type{MPlot}, m::Mesh, x) = (m, x, nothing)

# One parameter: Scalars or vector
function Makie.convert_arguments(::Type{MPlot}, m::Mesh, x::AbstractArray)

    veconly(u) = (m, LinearAlgebra.norm.(eachcol(u)), u)

    if ndims(x) == 1
        if length(x) == nnodes(m) || length(x) == nelements(m)
            return (m, x, nothing)
        end
        if length(x) == 2 * nnodes(m)
            return veconly(reshape(x, 2, :))
        end
    end

    if ndims(x) == 2 && size(x) == (2, nnodes(m))
        return veconly(x)
    end

    return (m, x, nothing)
end

function Makie.convert_arguments(::Type{MPlot}, m::Mesh, s::AbstractArray, u::AbstractArray)
    if ndims(u) == 1
        uu = reshape(u, 2, :)
    else
        uu = u
    end
    return (m, s, uu)
end

function Makie.plot!(plot::MPlot)
    _configure!(plot)
    plot.lineplotvisible[] && _plotlineplot(plot)
    plot.facesvisible[] && _plotfaces(plot)
    plot.faceplotfunction[] && _faceplotfunction(plot)
    plot.edgesvisible[] && _plotedges(plot, false)
    plot.featureedgesvisible[] && _plotedges(plot, true)
    plot.nodesvisible[] && _plotnodes(plot)

    if !isnan(plot.uscale[])
        Makie.text!(
            plot, 0, 1,
            text="×$(round(plot.uscale[], digits=4))",
            align=(:left, :top),
            offset=(4, -2),
            space=:relative
        )
    end

    return plot
end


# -------------------------------------------------------------------------------------------------
# Plot configuration
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
        if colorbar && !scene.plot.groupsforcolors[]
            p = _find_plot_with_colormap(scene.plot.plots)
            if !isnothing(p) &&
               !isnothing(Makie.extract_colormap_recursive(p)) &&
               !(p isa Makie.Lines)
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
# Plot logic
# -------------------------------------------------------------------------------------------------

function _configure!(plot::MPlot)

    # Mesh
    mesh = plot.mesh[]

    # Smooth
    if plot.smooth[]
        plot.scalars = _smooth(mesh, plot.scalars_raw[])
    else
        plot.scalars = plot.scalars_raw[]
    end

    # Shorthand
    havegroups = false
    havescalars = !isnothing(plot.scalars[])
    havevectors = !isnothing(plot.vectors[])

    # Visibility off
    if (nnodes(mesh)) == 0
        plot.nodesvisible = false
    end
    if nedges(mesh) == 0
        plot.edgesvisible = false
        plot.lineplotvisible = false
        plot.featureedgesvisible = false
    end
    if nfaces(mesh) == 0
        plot.facesvisible = false
    end

    # Warp node coordinates
    nwarp(node, _=1) = [coordinates(node)..., nodewarp[index(node)]]
    ndisp(node, scale=1) = coordinates(node) + scale * plot.vectors[][:, index(node)]
    nid(node, _=1) = coordinates(node)

    if _isdefined(plot, :nodewarp)
        nodewarp = plot.nodewarp[]
        if nodewarp isa Function
            plot.nodecoordinates = (n, _) -> plot.nodewarp[](n)
        elseif nodewarp isa AbstractVector
            plot.nodecoordinates = nwarp
        end
    elseif havevectors
        plot.nodecoordinates = ndisp
    else
        plot.nodecoordinates = nid
    end

    # Displacement scale
    if havevectors && !_isdefined(plot, :uscale)
        umax = maximum(norm.(eachcol(plot.vectors[])))
        size = mesh.geometry |> boundingbox |> diagonal
        plot.uscale = 0.1 * size / umax
    end

    # Settings for mesh of 1D elements
    if pdim(mesh) == 1
        plot.edgesvisible = true
        plot.facesvisible = false
        plot.featureedgesvisible = false
        _setifundefined(plot, :edgelinewidth, 3)
        _setifundefined(plot, :lineplotvisible, havescalars)

        if plot.groupsforcolors[] && hasgroups(mesh, d=1) && !havescalars
            havegroups = true
            _setifundefined(plot, :edgecolor, :groups)
        end
    end

    # # Settings for mesh of 2D elements
    if pdim(mesh) == 2
        _setifundefined(plot, :lineplotvisible, false, enforce=true)
        _setifundefined(plot, :edgesvisible, nfaces(mesh) <= 100)
        _setifundefined(plot, :edgelinewidth, DEFAULT_EDGE_LINEWIDTH_2D)
        _setifundefined(plot, :featureedgelinewidth, DEFAULT_FEATUREEDGE_LINEWIDTH_2D)

        if havescalars
            if plot.scalars[] isa Function || plot.scalars[] isa Symbol
                plot.nodesvisible = false
                plot.facesvisible = false
                plot.faceplotfunction = true
                _setifundefined(plot, :facecolor, DEFAULT_FACE_COLORCOMPONENT)
            end
        else
            if plot.groupsforcolors[]
                if hasgroups(mesh, d=1)
                    havegroups = true
                    _setifundefined(plot, :featureedgecolor, :groups)
                end
                if hasgroups(mesh, d=2)
                    havegroups = true
                    _setifundefined(plot, :facecolor, :groups)
                    _setifundefined(plot, :facecolormap, DEFAULT_FACE_COLORMAP_GROUPS)
                end
            end
        end
    end

    # Higher dimensions not implemented yet
    if pdim(mesh) >= 3
        @notimplemented
    end

    # Apply default values
    plot.groupsforcolors = havegroups
    _setifundefined(plot, :nodesvisible, nnodes(mesh) <= 50)
    _setifundefined(plot, :edgecolor, Makie.theme(:linecolor)[])
    _setifundefined(plot, :featureedgecolor, Makie.theme(:linecolor)[])
    _setifundefined(plot, :facecolor, DEFAULT_FACE_COLOR)
    _setifundefined(plot, :facecolormap, Makie.theme(:colormap)[])
    _setifundefined(plot, :uscale, NaN)
end


# -------------------------------------------------------------------------------------------------
# Plot parts
# -------------------------------------------------------------------------------------------------

_plotlineplot(plot::MPlot) =
    _plotlineplot(plot, plot.mesh[], plot.scalars[])

# Version for numbers
_plotlineplot(plot::MPlot, mesh::Mesh, values::AbstractVecOrMat{<:Real}) =
    _plotlineplot(plot, mesh, _tofunctions(mesh, values))

# Version for functions
function _plotlineplot(plot::MPlot, mesh::Mesh, functions::AbstractVector{<:FunctionRToR{IHat}})

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
    if !_isautomatic(plot, :lineplotfacescolor)
        cf = plot.lineplotfacescolor[]
    end

    # Plot faces
    if plot.lineplotfacesvisible[]
        Makie.mesh!(plot, [xf yf]', reshape(triangles, 3, :)',
            color=cf, colormap=plot.lineplotfacescolormap)
    end

    # Plot outlines
    if plot.lineplotoutlinesvisible[]
        Makie.lines!(plot, xe, ye, linewidth=plot.lineplotoutlineslinewidth,
            color=plot.lineplotoutlinescolor)
    end
end


function _plotfaces(plot::MPlot)
    # Mesh and color
    mesh = plot.mesh[]
    scalars = plot.scalars[]
    havescalars = !isnothing(scalars)
    color = havescalars ? scalars : plot.facecolor

    # Coordinates
    nodecoordinates = plot.nodecoordinates[]

    # Test if we have data and overall color
    haveData = color isa Vector

    # Helpers
    coords = nodecoordinates.(nodes(mesh), plot.uscale[])
    Nn = nnodes(mesh)
    Nf = nfaces(mesh)

    # Color by groups
    if plot.facecolor[] == :groups
        color = groupids(mesh, d=2, predefined=false)
        color = maximum(color) .- color
        haveData = true
    end

    # Collect
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
    Makie.mesh!(
        plot,
        tomatrix(coords), tomatrix(tf, ROWS),
        color=color, colormap=plot.facecolormap, colorrange=plot.colorrange
    )
end


function _plotedges(plot::MPlot, featureedges::Bool)
    mesh = plot.mesh[]
    nodecoordinates = plot.nodecoordinates[]

    if featureedges
        linewidth = plot.featureedgelinewidth[]
        color = plot.featureedgecolor[]
        indices = group(mesh, :boundaryedges)
        if color == :groups
            for n in groupnames(mesh, d=1, predefined=false)
                indices = indices ∪ group(mesh, n)
            end
        end
        uscale = 0
    else
        linewidth = plot.edgelinewidth[]
        color = plot.edgecolor[]
        indices = 1:nedges(mesh)
        uscale = plot.uscale[]
    end

    if color == :groups
        ids = groupids(mesh, d=1, predefined=true)
        color = reshape(repeat(ids[indices], 1, 3)', :)
    end

    points = [nodecoordinates.(nodes(edge(mesh, i)), uscale) for i = indices]

    Makie.lines!(
        plot, _collectlines(points)...,
        linewidth=linewidth, color=color,
        colormap=DEFAULT_EDGE_COLORMAP
    )
end

function _plotnodes(plot::MPlot)
    mesh = plot.mesh[]
    nodecoordinates = plot.nodecoordinates[]
    x = nodecoordinates.(nodes(mesh), plot.uscale[]) |> tomatrix

    Makie.scatter!(
        plot, x, color=plot.nodecolor, markersize=plot.nodesize
    )
end

function _faceplotfunction(plot::MPlot)

    # Edge attributes
    edgesvisible = plot.edgesvisible[]
    edgecolor = plot.edgecolor[]
    edgelinewidth = plot.edgelinewidth[]

    # Face plot
    facecolor = plot.facecolor[]

    # Faceplot
    fzscale = plot.faceplotzscale[]
    fnpoints = plot.faceplotnpoints[]
    fnmeshlines = plot.faceplotmesh[]
    fmeshcolor = plot.faceplotmeshcolor[]
    fmeshlinewidth = plot.faceplotmeshlinewidth[]

    # General
    colorrange = plot.colorrange[]
    colormap = plot.colormap[]

    # Parameters
    mesh = plot.mesh[]
    colors = plot.scalars[]

    # Function to plot
    if colors isa Symbol
        facefunction(face) = data(face, :post)(face, colors)
    else
        facefunction = colors
    end

    # Collect
    meshes, lpedges, lpmesh = _samplefaces(mesh, facefunction, fnpoints, fnmeshlines, fzscale)
    xf, tf = _mergemeshes(meshes)
    color = _getcolor(xf, facecolor, 1)
    xf[3, :] *= fzscale

    # Plot Mesh
    Makie.mesh!(
        plot, xf, tf, color=color,
        colormap=colormap, colorrange=colorrange
    )

    # Plot mesh lines on elements
    if fnmeshlines > 0
        Makie.lines!(plot, _collectlines(lpmesh)..., color=fmeshcolor, linewidth=fmeshlinewidth)
    end

    # Plot element edges
    if edgesvisible
        Makie.lines!(plot, _collectlines(lpedges)..., color=edgecolor, linewidth=edgelinewidth)
    end
end


# -------------------------------------------------------------------------------------------------
# Configure
# -------------------------------------------------------------------------------------------------

function _find_plot_with_colormap(plots)
    for p = plots
        !isnothing(Makie.extract_colormap_recursive(p)) && return p
    end
    return nothing
end

