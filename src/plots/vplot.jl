Makie.@recipe VPlot begin
    npoints = 10
    color = :black
    arrowsize = 7.5
    lengthscale = 0.25
    Makie.mixin_colormap_attributes()...
    Makie.mixin_generic_plot_attributes()...
end


function Makie.plot!(plot::VPlot)
    attributes = plot.attributes
    color = attributes.color[]
    np = attributes.npoints[]
    as = attributes.arrowsize[]
    ls = attributes.lengthscale[]
    cm = attributes.colormap[]

    for f ∈ plot.args[]
        d = domain(f)
        rx = DomainSets.component(d, 1)
        ry = DomainSets.component(d, 2)
        xs = range(rx, np)
        ys = range(ry, np)
        p = reshape([[x, y] for x in xs, y in ys], :)
        u = f.(p)
        if color == :norm
            color = norm.(u)
        end
        Makie.arrows2d!(
            plot, Makie.Point2f.(p), Makie.Vec2f.(u),
            # FIXME
            # linecolor=color,
            # arrowcolor=color, arrowsize=as, lengthscale=ls, colormap=cm
        )
    end
    return plot
end