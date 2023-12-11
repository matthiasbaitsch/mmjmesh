using Test
import CairoMakie

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Utilities
using MMJMesh.Plots

# 1D, horizontal
m = makemeshoninterval(0.0, 1.2, 10)
plot(m, -1.1 .+ 2.2 * rand(nnodes(m)))
@test true

# 1D, vertical
m = makemeshoninterval(π, 3π, 20, t -> [0; t])
plot(m, -1.1 .+ 3.2 * rand(nedges(m)))
@test true

# 2D
a = 3
m = makemeshonrectangle(9.0, 4.5, 2a, a)
plot(m, 4.1 * (rand(nnodes(m)) .- 0.25))
@test true
