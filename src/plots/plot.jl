"""
    PointConfiguration

How to plot points.
"""
mutable struct PointConfiguration
    visible::Bool
    size::Real
    color

    function PointConfiguration()
        return new(false, 7.0, :tomato)
    end
end

"""
    EdgeConfiguration

How to plot edges.
"""
mutable struct EdgeConfiguration
    visible::Bool
    outlineonly::Bool
    linewidth::Real
    color

    function EdgeConfiguration()
        return new(true, true, 1.0, :black)
    end
end

"""
    FaceConfiguration

How to plot faces.
"""
mutable struct FaceConfiguration
    visible::Bool
    color
    colormap

    function FaceConfiguration()
        return new(true, :lightsteelblue, :jet)
    end
end

"""
    PlotMeshConfiguration

How to plot the mesh.
"""
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
            true, # hidedecorations
            true  # colorbar
        )
    end
end

"""
    plot(m::Mesh[, pc::PlotMeshConfiguration=PlotMeshConfiguration()])

Plot mesh `m` with configuration `pc`.
"""
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
        Makie.hidespines!(ax)
        Makie.tightlimits!(ax)
    end

    # Return
    return f
end

function plotfaces(f::Makie.Figure, m::Mesh, cfg::FaceConfiguration, colorbar::Bool)

    # Test if we have data
    haveData = cfg.color isa Vector

    # Helpers
    coords = coordinates(m.geometry)
    nnodes = nentities(m.topology, 0)
    nfaces = nentities(m.topology, 2)

    # Collect triangles
    
    # No color or nodal color
    if !haveData || length(cfg.color) == nnodes
        t = Vector{Vector{Int}}()
        for l in links(m.topology, 2, 0)
            if length(l) == 3
                push!(t, l)
            elseif length(l) == 4
                push!(t, [l[1], l[2], l[3]])
                push!(t, [l[1], l[3], l[4]])
            end
        end
        x = coords
        t = mapreduce(permutedims, vcat, t)
        c = cfg.color
    # Element color
    elseif length(cfg.color) == nfaces
        cnt = 1;
        t = Vector{Vector{Int}}()
        c = Vector{Float64}()
        xc = Vector{Float64}()
        yc = Vector{Float64}()
        for (i, l) ∈ enumerate(links(m.topology, 2, 0))
            if length(l) == 3
                append!(c, cfg.color[i] * ones(3))
                append!(xc, coords[1, l[[1, 2, 3]]])
                append!(yc, coords[2, l[[1, 2, 3]]])
                push!(t, cnt:cnt+2)
                cnt += 3
            elseif length(l) == 4
                append!(c, cfg.color[i] * ones(6))
                append!(xc, coords[1, l[[1, 2, 3, 1, 3, 4]]])
                append!(yc, coords[2, l[[1, 2, 3, 1, 3, 4]]])
                push!(t, cnt:cnt+2)
                push!(t, cnt+3:cnt+5)
                cnt += 6
            else
                @notimplemented
            end
        end
        x = [xc yc]'
        t = mapreduce(permutedims, vcat, t)
    else
        @notimplemented
    end

    # Plot
    hm = Makie.mesh!(x, t, color=c, colormap=cfg.colormap)

    # Colorbar
    if haveData && colorbar
        Makie.Colorbar(f[1, 2], hm)
    end
end

function plotedges(m::Mesh, cfg::EdgeConfiguration)

    # Coordinates and links from edges to faces
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
