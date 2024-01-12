using Test
using MMJMesh.Topologies

# Basic
@testset "Basic" begin
    cl = ConnectivityList([[1, 3, 2], [1, 8, 4, 2], [1, 9, 3, 2], [3, 2]])
    push!(cl, [1, 2])
    push!(cl, [3, 4, 5])
    @test length(cl) == 6
    @test length(cl, 2) == 4
    @test maxlinkssize(cl) == 4
    @test cl[3] == [1, 9, 3, 2]
    @test cl[5] == [1, 2]
    @test cl[6] == [3, 4, 5]
    @test cl[3, 2] == 9
    @test [1, 3, 2] ∈ cl
    @test [3, 2] ∈ cl
    @test [1, 22] ∉ cl
end

# Transpose
@testset "Transpose" begin
    cl = ConnectivityList([[1, 3, 2], [1, 8, 4, 2], [1, 9, 3, 2], [3, 2]])
    clt = cl'
    @test length(clt) == 9
    @test clt[1] == [1, 2, 3]
    @test clt[2] == [1, 2, 3, 4]
    @test clt[7] == []
    @test clt[8] == [2]
    @test clt[9] == [3]
end

# Inverse
@testset "Inverse" begin
    cl = ConnectivityList([[1, 3, 2], [1, 8, 4, 2], [1, 9, 3, 2], [3, 2]])
    cli = inverse(cl)
    @test cli[Set([1, 3, 2])] == 1
    @test cli[Set([1, 9, 3, 2])] == 3
    @test cli[Set([3, 2])] == 4
end

# Comparison
@testset "Comparison" begin
    cl1 = ConnectivityList([[1, 2, 5, 4], [2, 3, 5], [3, 6, 5]])
    cl2 = ConnectivityList([[1, 2, 5, 4], [2, 3, 5], [3, 6, 5]])
    @test cl1 == cl2
end
