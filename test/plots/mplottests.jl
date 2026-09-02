# TODO scratch script for visually eyeballing mplot output during
# development -- has no @test assertions, so it isn't wired into the test
# suite (like test/plots/fplottests.jl). Either add assertions and convert
# to @testitem, or move it out of test/ alongside dev.jl.

using MMJMesh
using MMJMesh.Plots
using MMJMesh.Meshes

import CairoMakie

mplot(Mesh(0 .. 4, 20)) |> mconf()

