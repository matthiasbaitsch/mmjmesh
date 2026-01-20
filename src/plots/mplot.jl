# -------------------------------------------------------------------------------------------------
# Recipe
# -------------------------------------------------------------------------------------------------

Makie.@recipe MPlot (mesh, scalars) begin

    # Nodes
    nodesvisible = nothing
    nodecolor = :tomato
    nodesize = @inherit markersize
    nodewarp = nothing

    # Edges
    edgesvisible = nothing
    edgecolor = @inherit linecolor
    edgelinewidth = nothing

    # Featureedges
    featureedgesvisible = true
    featureedgecolor = nothing
    featureedgelinewidth = nothing

    # Faces
    facesvisible = true
    facecolor = nothing
    facecolormap = nothing

    # Lineplot
    lineplotvisible = nothing
    lineplotscale = 0.1
    lineplotoutlinesvisible = true
    lineplotoutlinescolor = :black
    lineplotoutlineslinewidth = @inherit linewidth
    lineplotfacesvisible = true
    lineplotfacescolor = nothing
    lineplotfacescolormap = @inherit colormap

    # Faceplot
    faceplotzscale = 0.0
    faceplotmesh = 15
    faceplotnpoints = 30
    faceplotmeshcolor = @inherit linecolor
    faceplotmeshlinewidth = 1.25

    # Defaults
    Makie.mixin_colormap_attributes()...
    Makie.mixin_generic_plot_attributes()...
end


# -------------------------------------------------------------------------------------------------
# Main plot function
# -------------------------------------------------------------------------------------------------

Makie.convert_arguments(::Type{MPlot}, m::Mesh) = (m, nothing)

function Makie.plot!(plot::MPlot)
    mesh = plot.mesh[]
    scalars = plot.scalars[]
    havescalars = !isnothing(scalars)

    # Helpers
    isundefined(key) = !haskey(plot, key) || isnothing(plot[key][])
    setifundefined(key, value) = isundefined(key) && (plot[key] = value)

    # Configure
    color = havescalars ? scalars : plot.facecolor
    colorbygroups = false
    lineplotvisible = false

    # Settings if not defined
    setifundefined(:nodesvisible, nnodes(mesh) <= 50)
    if (nnodes(mesh)) == 0
        plot.nodesvisible = false
    end
    if nedges(mesh) == 0
        plot.edgesvisible = false
        plot.featureedgesvisible = false
    end
    if nfaces(mesh) == 0
        plot.facesvisible = false
    end

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
        lineplotvisible = !isnothing(scalars)
        setifundefined(:edgelinewidth, 3)
    end

    # # Settings for mesh of 2D elements
    if pdim(mesh) == 2
        colorbygroups = !havescalars && isundefined(:facecolor) && hasgroups(mesh, d=2)
        setifundefined(:edgesvisible, nfaces(mesh) <= 100)
        setifundefined(:edgelinewidth, 0.75)
        setifundefined(:featureedgelinewidth, 2.25)

        if colorbygroups
            setifundefined(:facecolormap, :Pastel1_9)
        else
            if color isa Function || color isa Symbol
                setifundefined(:facecolor, 3)
            else
                setifundefined(:facecolor, :seashell2)
            end
            setifundefined(:facecolormap, Makie.theme(plot, :colormap)[])
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
        plotfeaturedges = plot.featureedgesvisible[] && nedges(mesh) > 0 && plot.faceplotzscale[] == 0

        plotfeaturedges && plotedges(plot, true)
        plot.facesvisible[] && plotfacefunctions(plot)
    else # Plot
        lineplotvisible && plotlineplot(plot)
        plot.facesvisible[] && plotfaces(plot, color)
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
        Makie.mesh!(plot, [xf yf]', reshape(triangles, 3, :)',
            color=cf, colormap=plot.lineplotfacescolormap)
    end

    # Plot outlines
    if plot.lineplotoutlinesvisible[]
        Makie.lines!(plot, xe, ye, linewidth=plot.lineplotoutlineslinewidth,
            color=plot.lineplotoutlinescolor)
    end
end


function plotfaces(plot::MPlot, color)

    # Mesh and color
    mesh = plot.mesh[]
    nodecoordinates = plot.nodecoordinates[]
    # color = length(plot) > 1 ? plot.scalars[] : plot.facecolor

    # Test if we have data and overall color
    haveData = color isa Vector
    haveColor = !isnothing(plot.facecolor[])

    # Helpers
    coords = nodecoordinates.(nodes(mesh))
    Nn = nnodes(mesh)
    Nf = nfaces(mesh)

    # Color by groups
    if !haveData && !haveColor && hasgroups(mesh, d=2)
        color = groupids(mesh, d=2, predefined=false)
        color = maximum(color) .- color
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
    Makie.mesh!(
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
        # Boundary edges
        indices = group(mesh, :boundaryedges)

        # Add edges in groups
        for n in groupnames(mesh, d=1, predefined=false)
            indices = indices ∪ group(mesh, n)
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
            if ngroups(mesh, d=1) > 0
                ids = groupids(mesh, d=1, predefined=true)
                lc = reshape(repeat(ids[indices], 1, 3)', :)
            else
                lc = :black
            end
        end
    else
        lw = plot.edgelinewidth[]
        lc = plot.edgecolor[]

        if lc == :groups
            @assert ngroups(mesh, d=1) > 0
            ids = groupids(mesh, d=1, predefined=true)
            lc = reshape(repeat(ids[indices], 1, 3)', :)
            plot.colorbygroups = true
        end
    end

    Makie.lines!(plot, _collectlines(points)..., linewidth=lw, color=lc, colormap=:tab10)
end

function plotnodes(plot::MPlot)
    mesh = plot.mesh[]
    nodecoordinates = plot.nodecoordinates[]
    x = nodecoordinates.(nodes(mesh)) |> tomatrix

    Makie.scatter!(
        plot, x, color=plot.nodecolor, markersize=plot.nodesize
    )
end

function plotfacefunctions(plot::MPlot)
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


