module MMJMeshRunTests

using Test
# using Aqua
using MMJMesh

@testset "MMJMesh" begin    
    # TODO Aqua.test_all(MMJMesh)
    @testset "Geometries" include("geometries/runtests.jl")
    @testset "Topologies" include("topologies/runtests.jl")    
    @testset "Associations" include("associations/00-runtests.jl")    
    @testset "Meshes" include("meshes/runtests.jl")
    @testset "Mathematics" include("mathematics/runtests.jl")    
    @testset "Utilities" include("utilities/00-runtests.jl")    
    @testset "Plots" include("plots/00-runtests.jl")    
end

end
