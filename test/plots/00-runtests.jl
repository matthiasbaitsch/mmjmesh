module PlotsTest

using Test
import CairoMakie

using MMJMesh
using MMJMesh.Plots
using MMJMesh.Meshes
using MMJMesh.Utilities

@testset "Basic" begin

    # 1D, horizontal
    m = makemeshoninterval(0.0, 1.2, 10)
    mplot(m, -1.1 .+ 2.2 * rand(nnodes(m)))
    @test true

    # 1D, vertical
    m = makemeshoninterval(π, 3π, 20, t -> [0; t])
    mplot(m, -1.1 .+ 3.2 * rand(nedges(m)))
    @test true

    # 2D
    a = 3
    m = makemeshonrectangle(9.0, 4.5, 2a, a)
    mplot(m, 4.1 * (rand(nnodes(m)) .- 0.25))
    @test true

end

@testset "Sample" include("sampleadaptivetests.jl")
@testset "Sample2" include("sampleadaptive2tests.jl")
@testset "Reference" include("referencetests.jl")

end
