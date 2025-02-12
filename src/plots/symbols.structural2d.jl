module Structural2D

using StaticArrays
using Base.Iterators
using Makie: Point2f, lines!, DataAspect, Figure, Axis, MoveTo, LineTo, ClosePath, BezierPath

using MMJMesh
using MMJMesh.MMJBase
using MMJMesh.Plots.Symbols
using MMJMesh.Plots.Symbols: PSymbol


const Constraint2D = StaticVector{3,Bool}


function makeconstraint2dsymbolmap()
    a = 1.0
    h1 = a * sqrt(3) / 2
    dd = 0.2 * h1
    h2 = h1 + dd

    fc(x) = x |> flatten |> collect

    hline(h, a=a) = [MoveTo(Point2f([-a / 2, h])), LineTo(Point2f([a / 2, h]))]

    triangle() = [
        MoveTo(Point2f(0, 0)),
        LineTo(Point2f(-a / 2, -h1)),
        LineTo(Point2f(a / 2, -h1)),
        ClosePath()
    ]

    function hatch(h)
        n = 8
        d = a / n
        p = [-a / 2, h]
        p1(i) = Point2f(p + i * [d, 0])
        p2(i) = Point2f(p + i * [d, 0] + [d, dd])
        return [[MoveTo(p1(i)), LineTo(p2(i))] for i = 0:n-1] |> fc
    end

    s010 = BezierPath([triangle(), hline(-h2)] |> fc)
    s110 = BezierPath([triangle(), hatch(-h2)] |> fc)
    s111 = BezierPath([hline(0), hatch(-dd)] |> fc)
    s011 = BezierPath([hline(0), hline(-dd), hatch(-2dd)] |> fc)

    symbols = Dict{Constraint2D,PSymbol}()
    e1 = [1, 0]
    e2 = [0, 1]

    symbols[SA[true, false, false]] = PSymbol{Constraint2D}(s010, [-e1, e1], [π / 2, 3π / 2])
    symbols[SA[false, true, false]] = PSymbol{Constraint2D}(s010, [e2, -e2], [0, π])
    symbols[SA[true, true, false]] = PSymbol{Constraint2D}(s110, [e2, -e1, -e2, e1], [0, π / 2, π, 3π / 2])
    symbols[SA[true, false, true]] = PSymbol{Constraint2D}(s011, [-e1, e1], [π / 2, 3π / 2])
    symbols[SA[false, true, true]] = PSymbol{Constraint2D}(s011, [e2, -e2], [0, π])
    symbols[SA[true, true, true]] = PSymbol{Constraint2D}(s111, [e1, e2, -e1, -e2], [3π / 2, 0, π / 2, π])

    return symbols
end


const constraint2dsymbolmap = makeconstraint2dsymbolmap()

const constraints2dsymbolstyle = Dict(
    :strokewidth => 1,
    :markersize => 15,
    :color => :transparent,
    :strokecolor => :black
)


MMJMesh.Plots.Symbols.symbolfor(c::Constraint2D) = constraint2dsymbolmap[c]
MMJMesh.Plots.Symbols.stylefor(::Type{Constraint2D}) = constraints2dsymbolstyle


function demo()
    function doit(sc, p, c, ds...)
        addsymbol!(sc, SVector{3,Bool}(collect(c) .== 't'), p, collect(ds))
        map(d -> lines!([p, p + 0.3 * d] |> tomatrix, color=:orange), ds)
        p[1] += 1
        (p[1] + 1) % 5 == 0 && (p[1] = 0; p[2] += 1)
        return p
    end

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=DataAspect())

    p = [0, 0]
    e1 = [1, 0]
    e2 = [0, 1]

    sc = SymbolCollection()
    p = doit(sc, p, "ttf", -e1, e2)
    p = doit(sc, p, "ftf", -e1, -e1 + e2)
    p = doit(sc, p, "tff", e1, e1 - e2)
    p = doit(sc, p, "tff", -e1, -e1 - e2)
    p = doit(sc, p, "ttf", -e1, e1, -e2)
    p = doit(sc, p, "ttf", -e1, e2, -e2)
    p = doit(sc, p, "ftt", e2)
    p = doit(sc, p, "tft", e1)
    p = doit(sc, p, "ttt", e1, e1 + e2)
    p = doit(sc, p, "ttt", -e2, -e1 - e2)
    p = doit(sc, p, "ttt", -e1)
    p = doit(sc, p, "ttt", -e2)

    stylefor(Constraint2D)[:color] = :lightgray
    stylefor(Constraint2D)[:strokecolor] = :green
    drawsymbols!(ax, sc)

    return fig
end

end