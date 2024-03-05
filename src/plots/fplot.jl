MakieCore.@recipe(FPlot, functions) do scene
    attr = MakieCore.Attributes(
        color=MakieCore.theme(scene, :linecolor),
        linewidth=MakieCore.theme(scene, :linewidth), 
        linestyle=nothing,
        cycle=[:color],
    )
    MakieCore.generic_plot_attributes!(attr,)
    MakieCore.colormap_attributes!(attr, MakieCore.theme(scene, :colormap))
    return attr
end


function MakieCore.plot!(plot::FPlot)
    # TODO: Plot options, sample options
    for a ∈ plot.args
        f = a[]
        a, b = endpoints(domain(f))
        xy = sampleadaptive(f, a, b)
        MakieCore.lines!(plot, xy)
        # MakieCore.lines!(plot, x, y; plot.attributes)
        # MakieCore.plot!(MakieCore.Lines, plot.attributes, x, y)
        # MakieCore.plot!(MakieCore.Lines, rand(10))
    end
    return plot
end

