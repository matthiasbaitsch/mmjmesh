mutable struct PointConfiguration
    visible::Bool
    size::Real
    color

    function PointConfiguration()
        return new(false, 7.0, :tomato)
    end
end

mutable struct EdgeConfiguration
    visible::Bool
    outlineonly::Bool
    width::Real
    color

    function EdgeConfiguration()
        return new(true, true, 1.0, :black)
    end
end

mutable struct PlotMeshConfiguration

    # Parts
    nodes::PointConfiguration
    edges::EdgeConfiguration

    # Appearance
    hidedecorations::Bool

    function PlotMeshConfiguration()
        return new(
            # Parts
            PointConfiguration(),
            EdgeConfiguration(),
            # Appearance
            false
        )
    end
end

function plotMesh(m::Mesh, pc::PlotMeshConfiguration=PlotMeshConfiguration())

    # Figure
    f = Figure()
    ax = Axis(f[1, 1], aspect=DataAspect())

    # Edges
    if pc.edges.visible
        cl12 = pc.edges.outlineonly ? links(m.topology, 1, 2) : ConnectivityList()

        for (i, l) in enumerate(links(m.topology, 1, 0))
            if !pc.edges.outlineonly || length(cl12, i) == 1
                println(i, " ", l)
            end
        end
    end

    # Nodes
    if pc.nodes.visible
        scatter!(
            coordinates(m.geometry),
            color=pc.nodes.color,
            markersize=pc.nodes.size
        )
    end

    # Edges

    # Appearance
    if pc.hidedecorations
        hidedecorations!(ax)
    end

    # Return
    return f
end
