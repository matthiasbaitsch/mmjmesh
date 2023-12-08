"""
    NodeStyle

How to plot points.
"""
mutable struct NodeStyle
    visible::Bool
    size::Real
    color

    function NodeStyle()
        return new(false, 7.0, :tomato)
    end
end

"""
    EdgeStyle

How to plot edges.
"""
mutable struct EdgeStyle
    visible::Bool
    outlineonly::Bool
    linewidth::Real
    color

    function EdgeStyle()
        return new(true, true, 1.0, :black)
    end
end

"""
    FaceStyle

How to plot faces.
"""
mutable struct FaceStyle
    visible::Bool
    color
    colormap

    function FaceStyle()
        return new(true, :lightsteelblue, :jet)
    end
end

"""
    LineplotStyle

How to plot values along edges.
"""
mutable struct LineplotStyle

    # Data
    values

    # Scaling factor
    scale::Float64

    # Parts
    outlines::EdgeStyle
    faces::FaceStyle


    function LineplotStyle()
        return new(
            # Data
            nothing,        # values
            0.1,            # scale
            # Parts
            EdgeStyle(),
            FaceStyle()
        )
    end
end

"""
    PlotStyle

How to plot the mesh.
"""
mutable struct PlotStyle

    # Parts
    nodes::Union{Nothing,NodeStyle}
    edges::Union{Nothing,EdgeStyle}
    faces::Union{Nothing,FaceStyle}
    lineplot::Union{Nothing,LineplotStyle}

    # Appearance
    hidedecorations::Bool
    colorbar::Bool
    # title::String

    function PlotStyle()
        return new(
            # Parts
            NodeStyle(),
            EdgeStyle(),
            FaceStyle(),
            LineplotStyle(),
            # Appearance
            true,  # hidedecorations
            true #,  # colorbar
            # ""
        )
    end
end

function PlotStyle(m::Mesh, values=nothing)
    ps = PlotStyle()

    # Only Edges
    if pdim(m) == 1
        ps.faces = nothing
        ps.nodes.visible = true
        ps.edges.linewidth = 2
        if !isnothing(values)
            ps.lineplot.values = values
        else
            ps.lineplot = nothing
        end

    # Faces
    elseif pdim(m) == 2
        ps.lineplot = nothing
        if !isnothing(values)
            ps.faces.color = values
        end
        ps.edges.outlineonly = nfaces(m) > 100
        
    # Edges and faces is all we got at the moment
    else
        @notimplemented
    end

    return ps
end

Base.show(io::IO, ::MIME"text/plain", ps::PlotStyle) = dump(io, ps)


