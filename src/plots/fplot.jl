"""
    fplot(f, [g, h, ...]; kwargs)

Plots the function `f` and optionally more functions.
"""
Makie.@recipe FPlot begin
    "Color"
    color = @inherit linecolor
    "Linewidth"
    linewidth = @inherit linewidth
    "Linestyle"
    linestyle = @inherit linestyle
    "Maximum number of recurstion steps"
    maxrecursion = 15
    "Maximum angle between segments in degrees"
    maxangle = 2.5
    "Initial number of sampling points"
    npoints = 5
    "Scaling factor in y-direction, affects evaluation of angle"
    yscale = 1.0
    "Inserts points where the polyline crosses the x-axis"
    ir = false
    "Connect lines at points of interest"
    connect_jumps = false
end


function Makie.plot!(plot::FPlot{<:Tuple{Vararg{<:AbstractMapping}}})
    attributes = plot.attributes

    # Process functions to plot
    for (i, f) = enumerate(plot.args[])

        # Shorthand to sample current function from a to b
        s1d(a, b) = MMJMesh.Plots.sample1d(f, a, b,
            maxrecursion=attributes.maxrecursion[],
            maxangle=attributes.maxangle[],
            npoints=attributes.npoints[],
            yscale=attributes.yscale[],
            ir=attributes.ir[]
        )

        # Points of interest and interval
        pts = pois(f)
        a, b = IntervalSets.endpoints(domain(f))

        # Collect x and y
        if isempty(pts) # No points of interest
            xy = s1d(a, b)
        else          # Handle points of intereste
            x = Float64[]
            y = Float64[]
            δ = 1e-12 * (b - a)

            # Shorthand to sample and to append to x and y
            s1da(a, b) = begin
                xy = s1d(a, b)
                append!(x, xy[1, :])
                append!(y, xy[2, :])

                if !attributes.connect_jumps[]
                    push!(x, NaN)
                    push!(y, NaN)
                end
            end

            # Sample first, inner and last intervals
            s1da(a, pts[1] - δ)
            for i = 1:length(pts)-1
                s1da(pts[i] + δ, pts[i+1] - δ)
            end
            s1da(pts[end] + δ, b)

            # Put together
            xy = stack([x, y], dims=1)
        end

        # Plot
        # Makie.lines!(plot, plot.attributes, xy; color=Makie.Cycled(i))
        Makie.lines!(plot, xy, color=Makie.Cycled(i))
    end

    return plot
end

fplot(functions::AbstractArray; kwargs...) = fplot(functions...; kwargs...)
