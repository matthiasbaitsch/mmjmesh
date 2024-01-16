using ReferenceTests
using Random
using CairoMakie
using MMJMesh
using MMJMesh.Plots
using MMJMesh.Groups
using MMJMesh.Meshes
using MMJMesh.Utilities
using MMJMesh.Gmsh

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

# Plot boundary nodes
a = 20
m = makemeshonrectangle(9.0, 4.5, 2a, a)
p = mplot(m)
x = coordinates(m)
bn = m.groups[:boundarynodes]
scatter!(p.axis, x[:, bn], color=:magenta)
@test_reference ref("m2d-008.png") p |> mconf()

# Face groups
a = 5
m = makemeshonrectangle(4, 2, 2a, a)
m.groups[:g1] = FaceGroup([1, 2, 3, 6, 22])
m.groups[:g2] = FaceGroup([5, 6, 7, 8, 22, 33])
m.groups[:g3] = FaceGroup([34])
@test_reference ref("m2d-009.png") mplot(m) |> mconf()
@test_reference ref("m2d-010.png") mplot(m, facecolor=:orange) |> mconf()

# Edge groups
a = 5
m = makemeshonrectangle(4, 2, 2a, a)
m.groups[:g1] = EdgeGroup(1:10)
m.groups[:g2] = EdgeGroup(8:16)
m.groups[:g3] = EdgeGroup(62:71)
@test_reference ref("m2d-011.png") mplot(m) |> mconf()
@test_reference ref("m2d-012.png") mplot(m, featureedgecolor=:orange) |> mconf()

# Gmesh meshes with groups
meshpath(m) = joinpath(@__DIR__(), "../../data/gmsh", m)
m = MMJMesh.Gmsh.Mesh(meshpath("advanced.msh"))
@test_reference ref("m2d-013.png") mplot(m) |> mconf()
m = MMJMesh.Gmsh.Mesh(meshpath("complex-g1.msh"))
@test_reference ref("m2d-014.png") mplot(m) |> mconf()
