Makie.@recipe(MPlot, mesh) do scene
    Makie.Attributes(
        
        # Nodes
        shownodes=nothing,
        nodecolor=nothing,
        nodesize=nothing,

        # Edges
        showedges=nothing,
        edgecolor=nothing,
        edgelinewidth=nothing,

        # Featureedges
        showfeatureedges=nothing,
        featureedgecolor=nothing,
        featureedgelinewidth=nothing,

        # Faces
        showfaces=nothing,
        facecolor=nothing,
        facecolormap=nothing,

        # Lineplot
        showlineplot=nothing,
        lineplotscale=nothing,
        showlineplotoutlines=nothing,
        lineplotoutlinescolor=nothing,
        lineplotoutlineslinewidth=nothing,
        showlineplotfaces=nothing,
        lineplotfacescolor=nothing,
        lineplotfacescolormap=nothing
    )
end

setoneattribute(s, key, value) = !isnothing(value[]) && setproperty!(s, key, value[])

# Transfer specified parameters to PlotStyle
function setattributes(s, plot::MPlot)

    if !isnothing(s.nodes)
        setoneattribute(s.nodes, :visible, plot.attributes.shownodes)
        setoneattribute(s.nodes, :color, plot.attributes.nodecolor)
        setoneattribute(s.nodes, :size, plot.attributes.nodesize)
    end

    if !isnothing(s.edges)
        setoneattribute(s.edges, :visible, plot.attributes.showedges)
        setoneattribute(s.edges, :color, plot.attributes.edgecolor)
        setoneattribute(s.edges, :linewidth, plot.attributes.edgelinewidth)
    end

    if !isnothing(s.featureedges)
        setoneattribute(s.featureedges, :visible, plot.attributes.showfeatureedges)
        setoneattribute(s.featureedges, :color, plot.attributes.featureedgecolor)
        setoneattribute(s.featureedges, :linewidth, plot.attributes.featureedgelinewidth)
    end

    if !isnothing(s.faces)
        setoneattribute(s.faces, :visible, plot.attributes.showfaces)
        setoneattribute(s.faces, :color, plot.attributes.facecolor)
        setoneattribute(s.faces, :colormap, plot.attributes.facecolormap)
    end

    if !isnothing(s.lineplot)
        setoneattribute(s.lineplot, :visible, plot.attributes.showlineplot)
        setoneattribute(s.lineplot, :scale, plot.attributes.lineplotscale)
        setoneattribute(s.lineplot.outlines, :visible, plot.attributes.showlineplotoutlines)
        setoneattribute(s.lineplot.outlines, :color, plot.attributes.lineplotoutlinescolor)
        setoneattribute(s.lineplot.outlines, :linewidth, plot.attributes.lineplotoutlineslinewidth)
        setoneattribute(s.lineplot.faces, :visible, plot.attributes.showlineplotfaces)
        setoneattribute(s.lineplot.faces, :color, plot.attributes.lineplotfacescolor)
        setoneattribute(s.lineplot.faces, :colormap, plot.attributes.lineplotfacescolormap)
    end
end

function Makie.plot!(plot::MPlot)
    # Use [] to get value from Observable

    # Mesh
    mesh = plot.mesh[]

    # Default plot style
    if length(plot) == 1                         # Only mesh
        s = PlotStyle(mesh)
    else                                         # Mesh and scalars
        scalars = plot[2][]
        s = PlotStyle(mesh, scalars)
    end

    # Set specified attributes
    setattributes(s, plot)

    # Do the plotting
    doplot(plot, mesh, s)

    # Done
    return plot
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
