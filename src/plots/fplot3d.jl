Makie.@recipe FPlot3D begin
    npoints = 30
    color = 3
    meshcolor = :black
    mesh = 5
    zscale = 1

    Makie.mixin_colormap_attributes()...
    Makie.mixin_generic_plot_attributes()...
end


function Makie.plot!(plot::FPlot3D)
    npoints = plot.npoints[]
    color = plot.color[]
    meshcolor = plot.meshcolor[]
    mesh = plot.mesh[]
    zscale = plot.zscale[]
    colorrange = plot.colorrange[]
    colormap = plot.colormap[]
    for f ∈ plot.args[]
        d = domain(f)
        !isfinite(d) && error("Domain is infinite: $d")
        fx, ft = sample2d(f, domain=d, npoints=npoints)
        fx[3, :] *= zscale
        lx = sample2dlines(f, domain=d, npoints=npoints, nmeshlines=mesh, zscale=zscale)
        Makie.mesh!(
            plot, fx, ft,
            color=_getcolor(fx, color, zscale), colorrange=colorrange, colormap=colormap
        )
        Makie.lines!(plot, _collectlines(lx)..., color=meshcolor)
    end
    return plot
end


function fplot3d(
    fs::AbstractArray{<:AbstractMapping}; colormap=Makie.theme(:colormap), fig=Makie.Figure(), zrange=nothing
)
    n = length(fs)
    cnt = 1
    ncol = Int(ceil(sqrt(n)))
    nrow = Int(ceil(n / ncol))
    for i = 1:ncol, j = 1:nrow
        if cnt <= n
            ax = Makie.Axis3(fig[i, j], protrusions=0)
            Makie.hidedecorations!(ax)
            fplot3d!(fs[cnt], colormap=colormap)
            !isnothing(zrange) && Makie.zlims!(ax, zrange)
            cnt += 1
        end
    end
    return fig
end
