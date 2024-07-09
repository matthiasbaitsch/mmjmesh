MakieCore.@recipe(VPlot, functions) do scene
    attr = Attributes(
        npoints=10,
        color=:black,
        arrowsize=7.5,
        lengthscale=0.25,
        colormap=MakieCore.theme(scene, :colormap)
    )
    MakieCore.generic_plot_attributes!(attr)
    return attr
end


function MakieCore.plot!(plot::VPlot)
    attributes = plot.attributes
    color = attributes.color[]
    np = attributes.npoints[]
    as = attributes.arrowsize[]
    ls = attributes.lengthscale[]
    cm = attributes.colormap[]

    for arg ∈ plot.args
        f = arg[]
        d = domain(f)
        rx = component(d, 1)
        ry = component(d, 2)
        xs = range(rx, np)
        ys = range(ry, np)
        p = reshape([[x, y] for x in xs, y in ys], :)
        u = f.(p)
        if color == :norm
            color = norm.(u)
        end
        arrows!(
            plot, Point2f.(p), Vec2f.(u), linecolor=color,
            arrowcolor=color, arrowsize=as, lengthscale=ls, colormap=cm
        )
    end
    return plot
end