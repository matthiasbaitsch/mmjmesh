module GeometryTests

using Test
using MMJMesh
using MMJMesh.Geometries
using MMJMesh.Geometries: coordinates

@testset "Empty constructor" begin
    g = Geometry(3)
    @test length(g, 0) == 0
    push!(g, [0.0 1.0; 0.0 0.0; 0.0 1.0])
    push!(g, [3.0, 0.0, 1.0])
    push!(g, [0.0, 1.0, 0.0])
    push!(g, [1.0, 1.0, 1.0])
    push!(g, [2.0, 1.0, 1.0])
    @test length(g, 0) == 6
    @test coordinates(g, 1) == [0.0, 0.0, 0.0]
    @test coordinates(g, [1, 3]) == [0.0 3.0; 0.0 0.0; 0.0 1.0]
end

@testset "Matrix constructor" begin
    g = Geometry([0 1 2; 0 0 9; 0 1 33])
    @test length(g, 0) == 3
    @test coordinates(g, 1) == [0.0, 0.0, 0.0]
    @test coordinates(g, [1, 3]) == [0.0 2.0; 0.0 9.0; 0.0 33.0]
end

@testset "Zeros constructor" begin
    g = Geometry(2, 5)
    g[0, 1] = Point(1, -1)
    g[0, 2] = [2 -6]
    g[0, 3] = [2, -2]
    g[0, 4] = 5:6
    @test g[0, 1] == Point(1, -1)
    @test g[0, 2] == Point(2, -6)
    @test g[0, 3] == Point(2, -2)
    @test g[0, 4] == Point(5, 6)
    @test g[0, 5] == Point(0, 0)
end

end