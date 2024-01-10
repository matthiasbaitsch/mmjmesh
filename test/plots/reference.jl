using ReferenceTests
using MMJMesh
using MMJMesh.Utilities
using MMJMesh.Plots
using MMJMesh.Meshes
using CairoMakie
using Random

Random.seed!(1234)
ref(f) = joinpath("../../data/references/plots", f)
update_theme!(colormap=:jet)

# -------------------------------------------------------------------------------------------------
# 1D meshes
# -------------------------------------------------------------------------------------------------
m = makemeshoninterval(0, 4, 60)
@test_reference ref("m1d-001.png") makemeshoninterval(0, 4, 20) |> mplot |> mconf()
@test_reference ref("m1d-002.png") m |> mplot |> mconf()
@test_reference ref("m1d-003.png") mplot(m, -1.1 .+ 2.6 * rand(nnodes(m))) |> mconf()
@test_reference ref("m1d-004.png") mplot(m, -1.1 .+ 2.2 * rand(nedges(m))) |> mconf()
@test_reference ref("m1d-005.png") mplot(m, -1.1 .+ 2.2 * rand(2, nedges(m))) |> mconf()
@test_reference ref("m1d-006.png") mplot(m, -1.1 .+ 2.2 * rand(2, nedges(m)),
    lineplotfacescolor=:gray50, lineplotoutlinescolor=:hotpink) |> mconf()


# -------------------------------------------------------------------------------------------------
# 2D meshes
# -------------------------------------------------------------------------------------------------

# Quadrilaterals
a = 80
m = makemeshonrectangle(9.0, 4.5, 2a, a)
@test_reference ref("m2d-001.png") mplot(m, edgesvisible=true, edgecolor=:hotpink) |> mconf()
@test_reference ref("m2d-002.png") mplot(m, 4.1 * (rand(nnodes(m)) .- 0.25)) |> mconf()
@test_reference ref("m2d-003.png") mplot(m, 4.1 * (rand(nfaces(m)) .- 0.25)) |> mconf()

# Triangles
a = 20
m = makemeshonrectangle(9.0, 4.5, 2a, a, TRIANGLE)
@test_reference ref("m2d-004.png") mplot(m, edgesvisible=true) |> mconf()
@test_reference ref("m2d-005.png") mplot(m, 4.1 * (rand(nnodes(m)) .- 0.25)) |> mconf()
@test_reference ref("m2d-006.png") mplot(m, 4.1 * (rand(nfaces(m)) .- 0.25)) |> mconf()

# Configuration
a = 10
m1 = makemeshonrectangle(4, 2, 2a, a)
@test_reference ref("m2d-007.png") mplot(m1, 3 * rand(nfaces(m1)),
    nodesvisible=true, nodecolor=:hotpink, nodesize=12,
    edgesvisible=true, edgecolor=:lightblue, edgelinewidth=3,
    featureedgecolor=:red, featureedgelinewidth=6,
    facecolormap=:bluesreds
) |> mconf()


# XXX
# a = 80
# m = makemeshonrectangle(9.0, 4.5, 2a, a)
# mplot(m)
# bn = m.groups[:boundarynodes]
# scatter!(coordinates(m)[:, bn], color=:orange)
