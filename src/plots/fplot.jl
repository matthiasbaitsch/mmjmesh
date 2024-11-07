MakieCore.@recipe(FPlot, functions) do scene
    attr = MakieCore.Attributes(
        color=MakieCore.theme(scene, :linecolor),
        linewidth=MakieCore.theme(scene, :linewidth),
        linestyle=nothing,
        maxrecursion=15,
        maxangle=2.5,
        npoints=5,
        yscale=1.0,
        ir=false,
        connect_jumps=false,
        cycle=[:color],
    )
    MakieCore.generic_plot_attributes!(attr,)
    MakieCore.colormap_attributes!(attr, MakieCore.theme(scene, :colormap))
    return attr
end


function MakieCore.plot!(plot::FPlot)

    att = plot.attributes

    # TODO: Plot options, sample options
    for a ∈ plot.args
        f = a[]
        a, b = endpoints(domain(f))
        p = pois(f)

        s1d(a, b) = sample1d(f, a, b,
            maxrecursion=att.maxrecursion[],
            maxangle=att.maxangle[],
            npoints=att.npoints[],
            yscale=att.yscale[],
            ir=att.ir[]
        )

        if isempty(p)
            xy = s1d(a, b)
        else
            δ = 1e-12 * (b - a)
            x = Float64[]
            y = Float64[]

            s1da(a, b) = begin
                xy = s1d(a, b)
                append!(x, xy[1, :])
                append!(y, xy[2, :])

                if !att.connect_jumps[]
                    push!(x, NaN)
                    push!(y, NaN)
                end
            end

            s1da(a, p[1] - δ)
            for i = 1:length(p)-1
                s1da(p[i] + δ, p[i+1] - δ)
            end
            s1da(p[end] + δ, b)

            xy = [x'; y']
        end


        MakieCore.lines!(plot, xy)
    end
    return plot
end


fplot(functions::AbstractArray; kwargs...) = fplot(functions...; kwargs)
