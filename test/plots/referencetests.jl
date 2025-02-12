using Random
using IntervalSets
using DomainSets: ×
using ReferenceTests
import CairoMakie as cm
import GLMakie as gm

using MMJMesh
using MMJMesh.Gmsh
using MMJMesh.Plots
using MMJMesh.Meshes
using MMJMesh.Utilities
using MMJMesh.Mathematics
using MMJMesh.Plots.Symbols.Structural2D

# Make sure to use CairoMakie
cm.activate!()


# -------------------------------------------------------------------------------------------------
# Set up
# -------------------------------------------------------------------------------------------------
Random.seed!(1234)
cm.update_theme!(colormap=:jet)
cm.update_theme!(px_per_unit=2)


# -------------------------------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------------------------------
ref(f) = joinpath(@__DIR__(), "../../data/references/plots", f)
meshpath(m) = joinpath(@__DIR__(), "../../data/gmsh", m)


# -------------------------------------------------------------------------------------------------
# 1D meshes
# -------------------------------------------------------------------------------------------------
m = Mesh(0 .. 4, 60)
@test_reference ref("m1d-001.png") Mesh(0 .. 4, 20) |> mplot |> mconf()
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
m = Mesh((0 .. 9.0) × (0 .. 4.5), 2a, a)
@test_reference ref("m2d-001.png") mplot(m, edgesvisible=true, edgecolor=:hotpink) |> mconf()
@test_reference ref("m2d-002.png") mplot(m, 4.1 * (rand(nnodes(m)) .- 0.25)) |> mconf()
@test_reference ref("m2d-003.png") mplot(m, 4.1 * (rand(nfaces(m)) .- 0.25)) |> mconf()

# Triangles
a = 20
m = Mesh((0 .. 9.0) × (0 .. 4.5), 2a, a, TRIANGLE)
@test_reference ref("m2d-004.png") mplot(m, edgesvisible=true) |> mconf()
@test_reference ref("m2d-005.png") mplot(m, 4.1 * (rand(nnodes(m)) .- 0.25)) |> mconf()
@test_reference ref("m2d-006.png") mplot(m, 4.1 * (rand(nfaces(m)) .- 0.25)) |> mconf()

# Configuration
a = 10
m1 = Mesh((0 .. 4) × (0 .. 2), 2a, a)
@test_reference ref("m2d-007.png") mplot(m1, 3 * rand(nfaces(m1)),
    nodesvisible=true, nodecolor=:hotpink, nodesize=12,
    edgesvisible=true, edgecolor=:lightblue, edgelinewidth=3,
    featureedgecolor=:red, featureedgelinewidth=6,
    facecolormap=:bluesreds
) |> mconf()

# Plot boundary nodes
a = 20
m = Mesh((0 .. 9.0) × (0 .. 4.5), 2a, a)
p = mplot(m)
x = coordinates(m)
bn = m.groups[:boundarynodes]
cm.scatter!(p.axis, x[:, bn], color=:magenta)
@test_reference ref("m2d-008.png") p |> mconf()

# Face groups
a = 5
m = Mesh((0 .. 9.0) × (0 .. 4.5), 2a, a)
definegroup!(m, 2, :g1, [1, 2, 3, 6, 22])
definegroup!(m, 2, :g2, [5, 6, 7, 8, 22, 33])
definegroup!(m, 2, :g3, [34])


@test_reference ref("m2d-009.png") mplot(m) |> mconf()
@test_reference ref("m2d-010.png") mplot(m, facecolor=:orange) |> mconf()

# Edge groups
a = 5
m = Mesh((0 .. 4) × (0 .. 2), 2a, a)
definegroup!(m, 1, :g1, 1:10)
definegroup!(m, 1, :g2, 8:16)
definegroup!(m, 1, :g3, 62:71)

@test_reference ref("m2d-011.png") mplot(m) |> mconf()
@test_reference ref("m2d-012.png") mplot(m, featureedgecolor=:orange) |> mconf()

# Gmesh meshes with groups
m = Mesh(meshpath("advanced.msh"))
@test_reference ref("m2d-013.png") mplot(m) |> mconf()
m = Mesh(meshpath("complex-g1.msh"))
@test_reference ref("m2d-014.png") mplot(m) |> mconf()

# Gmesh with nodes on edge group
m = Mesh(meshpath("multi_lambda.msh"))
p = mplot(m, edgesvisible=true) |> mconf()
cm.scatter!(p.axis, coordinates(m)[:, m.groups[:ΓD0]])
@test_reference ref("m2d-015.png") p |> mconf()
@test_reference ref("m2d-016.png") mplot(m, rand(nnodes(m))) |> mconf()

# Here are 3D plots
gm.activate!()

# Warp in z-direction
a = 4
m = Mesh((0 .. 4) × (0 .. 2), 2a, a)
function warp(node)
    x = coordinates(node)
    x3 = 0.1 * sin(0.25 * pi * (x[1] - 2)) * sin(0.5 * pi * (x[2] - 1))
    return [x..., x3]
end

f = gm.Figure() # Warp by function
gm.Axis3(f[1, 1], aspect=:data)
mplot!(m, rand(nfaces(m)), nodewarp=warp)
@test_reference ref("m2d-017.png") f

f = gm.Figure() # Warp by nodal values
gm.Axis3(f[1, 1], aspect=:data)
mplot!(m, rand(nfaces(m)), nodewarp=0.5 * rand(nnodes(m)))
@test_reference ref("m2d-018.png") f

# Plot function on face
m = Mesh((0 .. 8) × (0 .. 4), 4, 2)
w(face) = x -> index(face) * (1 - x[1]^2) * (1 - x[2]^2)

f = gm.Figure() # Warp by one function on face
gm.Axis3(f[1, 1], aspect=:data)
mplot!(m, w, faceplotzscale=0.2, faceplotmesh=2)
@test_reference ref("m2d-019.png") f

function results(face, name) # Warp using postprocessing function
    if name == :w
        return x -> index(face) * (1 - x[1]^2) * (1 - x[2]^2)
    elseif name == :sigma
        s = Polynomial([0, π])
        return ProductFunction(Sin() ∘ s, Cos() ∘ s)
    end
end
setdata!(m, :post, results)

f = gm.Figure()
gm.Axis3(f[1, 1], aspect=:data)
mplot!(m, :w, faceplotzscale=0.2, faceplotmesh=2)
@test_reference ref("m2d-020.png") f

f = gm.Figure()
gm.Axis3(f[1, 1], aspect=:data)
mplot!(m, :sigma, faceplotzscale=0.2, faceplotmesh=2)
@test_reference ref("m2d-21.png") f

f = gm.Figure()
gm.Axis3(f[1, 1], aspect=:data)
mplot!(m, :sigma, faceplotzscale=0.2, faceplotmesh=2, facecolor=:tomato)
@test_reference ref("m2d-22.png") f


# -------------------------------------------------------------------------------------------------
# Symbols
# -------------------------------------------------------------------------------------------------

cm.activate!()
f = MMJMesh.Plots.Symbols.Structural2D.demo()
@test_reference ref("sym-01.png") f


# -------------------------------------------------------------------------------------------------
# Other
# -------------------------------------------------------------------------------------------------

# Plot function R2 -> R
ff = MPolynomial([0 2 0; 0 0 2], [1, -1, -1], QHat)
f = gm.Figure()
gm.Axis3(f[1, 1], aspect=:data)
fplot3d!(ff, zscale=0.5)
@test_reference ref("plot-01.png") f

f = gm.Figure()
gm.Axis3(f[1, 1], aspect=:data)
fplot3d!(ff, zscale=0.5, mesh=nothing)
@test_reference ref("plot-02.png") f
