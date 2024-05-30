import MMJMesh.Mathematics: domain
import MakieCore: mesh!, lines!, Attributes


function _getcolor(x::Matrix, color)
    if typeof(color) == Int && 1 <= color <= 3
        return x[color, :]
    end
    return color
end


MakieCore.@recipe(FPlot3D, functions) do scene
    attr = Attributes(
        npoints=30,
        color=3,
        meshcolor=:black,
        mesh=5
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

    for arg ∈ plot.args
        f = arg[]
        d = domain(f)
        x, t = sample2d(f, domain=d, npoints=npoints)
        l1, l2, l3 = sample2dlines(d, f, npoints, mesh)
        mesh!(plot, x, t, color=_getcolor(x, color))
        lines!(plot, l1, l2, l3, color=meshcolor)
    end
    return plot
end


function fplot3d(fs::AbstractArray{<:AbstractMapping})
    n = length(fs)
    fig = Makie.Figure()
    cnt = 1
    ncol = Int(ceil(sqrt(n)))
    nrow = Int(ceil(n / ncol))
    for i = 1:ncol, j = 1:nrow
        if cnt <= n
            Makie.hidedecorations!(Makie.Axis3(fig[i, j], protrusions=0))
            fplot3d!(fs[cnt])
            cnt += 1
        end
    end
    return fig
end
