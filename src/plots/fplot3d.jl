import MMJMesh.Mathematics: domain
import MakieCore: mesh!, lines!, Attributes


MakieCore.@recipe(FPlot3D, functions) do scene
    attr = Attributes(
        npoints=30,
        color=3,
        meshcolor=:black,
        mesh=5,
        gmap=x -> x,
        colorrange=MakieCore.automatic,
        colormap=MakieCore.theme(scene, :colormap),
        zscale=1
    )
    MakieCore.generic_plot_attributes!(attr)
    MakieCore.colormap_attributes!(attr, MakieCore.theme(scene, :colormap))
    return attr
end


function MakieCore.plot!(plot::FPlot3D)
    attributes = plot.attributes
    mesh = attributes.mesh[]
    color = attributes.color[]
    meshcolor = attributes.meshcolor[]
    npoints = attributes.npoints[]
    gmap = attributes.gmap[]
    colorrange = attributes.colorrange[]
    colormap = attributes.colormap[]
    zscale = attributes.zscale[]
    for arg ∈ plot.args
        f = arg[]
        d = domain(f)
        !isfinite(d) && error("Domain is infinite: $d")
        fx, ft = sample2d(f, domain=d, npoints=npoints, gmap=gmap)
        fx[3, :] *= zscale
        lx = sample2dlines(f, domain=d, npoints=npoints, nmeshlines=mesh, gmap=gmap, zscale=zscale)
        mesh!(
            plot, fx, ft,
            color=_getcolor(fx, color, zscale), colorrange=colorrange, colormap=colormap
        )
        lines!(plot, _collectlines(lx)..., color=meshcolor)
    end
    return plot
end


function fplot3d(
    fs::AbstractArray{<:AbstractMapping}; colormap=MakieCore.theme(:colormap), fig=Makie.Figure()
)
    n = length(fs)
    cnt = 1
    ncol = Int(ceil(sqrt(n)))
    nrow = Int(ceil(n / ncol))
    for i = 1:ncol, j = 1:nrow
        if cnt <= n
            Makie.hidedecorations!(Makie.Axis3(fig[i, j], protrusions=0))
            fplot3d!(fs[cnt], colormap=colormap)
            cnt += 1
        end
    end
    return fig
end
