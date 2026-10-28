Makie.@recipe FEPlot (element,) begin

end

function plotit!(plot, K::IntervalSets.AbstractInterval)
    p1 = IntervalSets.leftendpoint(K)
    p2 = IntervalSets.rightendpoint(K)
    Makie.lines!(plot, [p1, p2], [0, 0], color=:black, linewidth=3)
end

function plotit!(plot, K::Union{DomainSets.Rectangle,DomainSets.EuclideanUnitSimplex{2}})
    corners = points(K, :corners)
    Makie.poly!(plot, corners, color=:seashell, strokecolor=:black, strokewidth=3)
end

function Makie.plot!(plot::FEPlot)
    element = plot.element[]

    plotit!(plot, element.K)

    # The linear forms
    previous = nothing
    pValueAt = []
    pDerivativeAt = []
    pMixedDerivativeAt = []

    # Helpers
    p2(x::Real) = Makie.Point2f(x, 0)
    p2(x::AbstractArray) = Makie.Point2f(x)

    handle(lf::ValueAtLF) = push!(pValueAt, p2(lf.x))
    handle(lf::DerivativeAtLF) = push!(pDerivativeAt, p2(lf.x))

    function handle(lf::PDerivativeAtLF)
        if lf.n == [1, 1]
            push!(pMixedDerivativeAt, p2(lf.x))
        elseif lf.n == [0, 1]
            @assert previous.n == [1, 0]
            push!(pDerivativeAt, p2(lf.x))
        end
    end

    # Process
    for lf = element.N
        handle(lf)
        previous = lf
    end

    # Plot
    if !isempty(pValueAt)
        Makie.scatter!(plot, pValueAt, color=:black, markersize=13)
    end

    if !isempty(pDerivativeAt)
        Makie.scatter!(
            plot, pDerivativeAt,
            color=:transparent, strokewidth=1.5, strokecolor=:black, markersize=20
        )
    end

    if !isempty(pMixedDerivativeAt)
        Makie.arrows2d!(
            plot, 
            pMixedDerivativeAt,
            repeat([Makie.Vec2f(1, 1)], length(pMixedDerivativeAt)), 
            lengthscale=0.2
        )
    end

    return plot
end


"""
    _bounds(K) -> (xmin, xmax, ymin, ymax)

Bounding box of domain `K` for plotting. Intervals are drawn on the x-axis and get a square box.
"""
function _bounds(K::IntervalSets.AbstractInterval)
    a, b = IntervalSets.endpoints(K)
    return a, b, -(b - a) / 2, (b - a) / 2
end

function _bounds(K::Union{DomainSets.Rectangle,DomainSets.EuclideanUnitSimplex{2}})
    xs = [p[1] for p in points(K, :corners)]
    ys = [p[2] for p in points(K, :corners)]
    return minimum(xs), maximum(xs), minimum(ys), maximum(ys)
end

function feconf()
    return scene -> begin
        xmin, xmax, ymin, ymax = _bounds(scene.plot.element[].K)
        g = 0.15 * max(xmax - xmin, ymax - ymin)
        ax = scene.axis
        ax.aspect = Makie.DataAspect()
        Makie.xlims!(ax, xmin - g, xmax + g)
        Makie.ylims!(ax, ymin - g, ymax + g)
        Makie.hidedecorations!(ax)
        Makie.hidespines!(ax)
        return scene
    end
end
