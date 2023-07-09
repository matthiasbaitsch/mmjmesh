module TopologyTests

using Test
using MMJMesh

t = Topology(2, 6)
addlinks!(t, 2, 0, [[1, 2, 5, 4], [5, 2, 3], [6, 5, 3]])

@testset "Basic" begin
    @test nentities(t, 0) == 6
    @test entity(t, 0, 2) == 2
    @test nlinks(t, 1, 0) == 0
    @test nlinks(t, 2, 0) == 3
    @test links(t, 2, 0) == ConnectivityList([[1, 2, 5, 4], [5, 2, 3], [6, 5, 3]])
    @test link(t, 2, 0, 1) == [1, 2, 5, 4]
end

@testset "Node -> Face" begin
    cl = links(t, 0, 2)
    @test length(cl) == 6
    @test cl[1] == [1]
    @test cl[2] == [1, 2]
    @test cl[3] == [2, 3]
    @test cl[4] == [1]
    @test cl[5] == [1, 2, 3]
    @test cl[6] == [3]
end

@testset "Edge -> Node" begin
    cl = links(t, 1, 0)
    @test length(cl) == 8
    @test cl[1] == [1, 2]
    @test cl[2] == [2, 5]
    @test cl[3] == [4, 5]
    @test cl[4] == [1, 4]
    @test cl[5] == [2, 3]
    @test cl[6] == [5, 3]
    @test cl[7] == [6, 5]
    @test cl[8] == [6, 3]
end

@testset "Face -> Edge" begin
    cl = links(t, 2, 1)
    @test length(cl) == 3
    @test cl[1] == [1, 2, 3, 4]
    @test cl[2] == [2, 5, 6]
    @test cl[3] == [7, 6, 8]
end

@testset "Node -> Node" begin
    cl = links(t, 0, 0)
    @test length(cl) == 6
end

@testset "Node -> Node" begin
    cl = links(t, 0, 0)
    @test length(cl) == 6
    @test length(cl, 1) == 0
    @test length(cl, 2) == 0
    @test length(cl, 6) == 0    
end

@testset "Edge -> Edge" begin
    cl = links(t, 1, 1)
    @test length(cl) == 8
    @test cl[1] == [4, 2, 5]
    @test cl[2] == [1, 5, 3, 6, 7]
    @test cl[3] == [4, 2, 6, 7]
    @test cl[4] == [1, 3]
    @test cl[5] == [1, 2, 6, 8]
    @test cl[6] == [5, 8, 2, 3, 7]
    @test cl[7] == [2, 3, 6, 8]
    @test cl[8] == [5, 6, 7]
end

@testset "Face -> Face" begin
    cl = links(t, 2, 2)
    @test length(cl) == 3
    @test cl[1] == [2]
    @test cl[2] == [1, 3]
    @test cl[3] == [2]
end

end