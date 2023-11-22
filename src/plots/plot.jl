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
    linewidth::Real
    color

    function EdgeConfiguration()
        return new(true, true, 1.0, :black)
    end
end

mutable struct FaceConfiguration
    visible::Bool
    color
    colormap

    function FaceConfiguration()
        return new(true, :lightsteelblue, :jet)
    end
end

mutable struct PlotMeshConfiguration

    # Parts
    nodes::PointConfiguration
    edges::EdgeConfiguration
    faces::FaceConfiguration

    # Appearance
    hidedecorations::Bool
    colorbar::Bool

    function PlotMeshConfiguration()
        return new(
            # Parts
            PointConfiguration(),
            EdgeConfiguration(),
            FaceConfiguration(),
            # Appearance
            false,
            true
        )
    end
end

function plot(m::Mesh, pc::PlotMeshConfiguration=PlotMeshConfiguration())

    # Figure
    f = Makie.Figure()
    ax = Makie.Axis(f[1, 1], aspect=Makie.DataAspect())

    # Faces
    if pc.faces.visible
        plotfaces(f, m, pc.faces, pc.colorbar)
    end

    # Edges
    if pc.edges.visible
        plotedges(m, pc.edges)
    end

    # Nodes
    if pc.nodes.visible
        Makie.scatter!(
            coordinates(m.geometry),
            color=pc.nodes.color,
            markersize=pc.nodes.size
        )
    end

    # Appearance
    if pc.hidedecorations
        Makie.hidedecorations!(ax)
        Makie.tightlimits!(ax)
    end

    # Return
    return f
end

function plotfaces(f::Makie.Figure, m::Mesh, cfg::FaceConfiguration, colorbar::Bool)
    t = Vector{Vector{Int}}()
    for l in links(m.topology, 2, 0)
        if length(l) == 3
            push!(t, l)
        elseif length(l) == 4
            push!(t, [l[1], l[2], l[3]])
            push!(t, [l[1], l[3], l[4]])
        end
    end
    hm = Makie.mesh!(
        coordinates(m.geometry),
        mapreduce(permutedims, vcat, t),
        color=cfg.color,
        colormap=cfg.colormap
    )
    if colorbar
        Makie.Colorbar(f[1, 2], hm)
    end
end

function plotedges(m::Mesh, cfg::EdgeConfiguration)
    coords = coordinates(m.geometry)
    l12 = cfg.outlineonly ? links(m.topology, 1, 2) : ConnectivityList()

    # Collect coordinates - currently only 2D
    xx = Vector{Float64}()
    yy = Vector{Float64}()
    for (i, l) in enumerate(links(m.topology, 1, 0))
        if !cfg.outlineonly || length(l12, i) == 1
            x1 = coords[:, l[1]]
            x2 = coords[:, l[2]]
            push!(xx, x1[1], x2[1], NaN)
            push!(yy, x1[2], x2[2], NaN)
        end
    end

    # Plot
    Makie.lines!(xx, yy, linewidth=cfg.linewidth, color=cfg.color)
end
