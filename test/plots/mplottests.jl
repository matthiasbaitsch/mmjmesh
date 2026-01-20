using MMJMesh
using MMJMesh.Plots
using MMJMesh.Meshes

import CairoMakie

mplot(Mesh(0 .. 4, 20)) |> mconf()

