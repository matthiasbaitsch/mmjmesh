@kwdef mutable struct PointStyle
    visible::Bool = false
    size::Real = 7.0
    color = :tomato
end

@kwdef mutable struct LineStyle
    visible::Bool = true
    linewidth::Real = 1.0
    color = :gray30
end

@kwdef mutable struct MeshStyle
    visible::Bool = true
    color = :seashell2                         # Color or scalars
    colormap = :jet                             # Unused unless color is scalars
end

@kwdef mutable struct LineplotStyle
    visible = false
    values = nothing
    scale::Float64 = 0.1
    outlines::LineStyle = LineStyle()
    faces::MeshStyle = MeshStyle(color=nothing)
end

@kwdef mutable struct PlotStyle
    nodes::Union{Nothing,PointStyle} = PointStyle()
    edges::Union{Nothing,LineStyle} = LineStyle()
    featureedges::Union{Nothing,LineStyle} = LineStyle(color=:gray30)
    faces::Union{Nothing,MeshStyle} = MeshStyle()
    lineplot::Union{Nothing,LineplotStyle} = LineplotStyle()
end

function PlotStyle(m::Mesh, scalars=nothing)
    ps = PlotStyle()

    if pdim(m) == 1                              # 1D mesh
        ps.nodes.visible = true
        ps.edges.linewidth = 2.0
        ps.faces = nothing
        ps.featureedges = nothing
        if !isnothing(scalars)
            ps.lineplot.visible = true
            ps.lineplot.values = scalars
        else
            ps.lineplot = nothing
        end
    elseif pdim(m) == 2                          # 2D mesh
        ps.lineplot = nothing
        ps.edges.color = :gray60
        ps.featureedges.linewidth = 1.5
        if !isnothing(scalars)
            ps.faces.color = scalars
        end
        ps.edges.visible = nfaces(m) <= 100
    else
        @notimplemented
    end

    return ps
end

Base.show(io::IO, ::MIME"text/plain", ps::PlotStyle) = dump(io, ps)


