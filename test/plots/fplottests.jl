# TODO scratch script for visually eyeballing fplot output during
# development -- has no @test assertions, so it isn't wired into the test
# suite (like test/plots/mplottests.jl). Either add assertions and convert
# to @testitem, or move it out of test/ alongside dev.jl.

using MMJMesh
using MMJMesh.Plots
using MMJMesh.Mathematics

import CairoMakie

x = parameter(0..3π)

fplot(sin(x))

