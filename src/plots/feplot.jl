
MakieCore.@recipe(FEPlot, element) do scene
    attr = MakieCore.Attributes()
    return attr
end

function plotit!(plot, K::IntervalSets.AbstractInterval)
    p1 = IntervalSets.leftendpoint(K)
    p2 = IntervalSets.rightendpoint(K)
    MakieCore.lines!(plot, [p1, p2], [0, 0], color=:black, linewidth=3)
end

function plotit!(plot, K::DomainSets.Rectangle)
    corners = points(K, :corners)
    MakieCore.poly!(plot, corners, color=:seashell, strokecolor=:black, strokewidth=3)
end

function MakieCore.plot!(plot::FEPlot)
    element = plot.element[]

    plotit!(plot, element.K)

    # The linear forms
    previous = nothing
    pValueAt = []
    pDerivativeAt = []
    pMixedDerivativeAt = []

    # Helpers
    p2(x::Real) = MakieCore.Point2f(x, 0)
    p2(x::AbstractArray) = MakieCore.Point2f(x)

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
        MakieCore.scatter!(plot, pValueAt, color=:black, markersize=13)
    end

    if !isempty(pDerivativeAt)
        MakieCore.scatter!(
            plot, pDerivativeAt,
            color=:transparent, strokewidth=1.5, strokecolor=:black, markersize=20
        )
    end

    if !isempty(pMixedDerivativeAt)
        MakieCore.arrows2d!(
            plot, 
            pMixedDerivativeAt,
            repeat([MakieCore.Vec2f(1, 1)], length(pMixedDerivativeAt)), 
            lengthscale=0.2
        )
    end

    return plot
end


function feconf()
    return scene -> begin
        g = 0.3
        ax = scene.axis
        ax.aspect = Makie.DataAspect()
        Makie.xlims!(ax, -1 - g, 1 + g)
        Makie.ylims!(ax, -1 - g, 1 + g)
        Makie.hidedecorations!(ax)
        Makie.hidespines!(ax)
        return scene
    end
end