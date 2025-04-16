using MMJMesh
using MMJMesh.Gmsh
using MMJMesh.Plots
using MMJMesh.Meshes
using MMJMesh.Mathematics

using LinearAlgebra

import GLMakie
using GLMakie: Figure, Axis3

GLMakie.activate!()


dofindices(idxs, nf) = vcat([nf * (i - 1) .+ (1:nf) for i = idxs]...)

function make_post(uhat)
    N4 = makeelement(:lagrange, QHat, k=1) |> nodalbasis
    return (face, name) -> begin
        nidxs = nodeindices(face)
        ndofs = dofindices(nidxs, 3)
        ue = uhat[ndofs]
        name == :w && return ue[1:3:end] ⋅ N4
    end
end


m = Mesh("data/gmsh/paganetty-plate.geo")

setdata!(m, :post, make_post(rand(3 * nnodes(m))))
f = Figure()
Axis3(f[1, 1], aspect=:data)
@time mplot!(m, :w, faceplotzscale=1, faceplotnpoints=1)
f


