module MMJMeshRunTests

using Test

@testset "MMJMesh" begin    
    @testset "Geometries" include("geometries/runtests.jl")
    @testset "Topologies" include("topologies/runtests.jl")    
    @testset "Mathematics" include("mathematics/runtests.jl")    
end

end
