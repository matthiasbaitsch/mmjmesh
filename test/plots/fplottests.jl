using MMJMesh
using MMJMesh.Plots
using MMJMesh.Mathematics

import CairoMakie

x = parameter(0..3π)

fplot(sin(x))

