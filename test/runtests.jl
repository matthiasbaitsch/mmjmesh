module MMJMeshRunTests

using Test

@testset "MMJMesh" begin    
    @testset "Geometries" include("geometries/runtests.jl")
    @testset "Topologies" include("topologies/runtests.jl")    
end

end
