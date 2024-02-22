module MMJMeshRunTests

using Test
using MMJMesh

@testset "MMJMesh" begin    
    @testset "Geometries" include("geometries/00-runtests.jl")
    @testset "Topologies" include("topologies/00-runtests.jl")    
    @testset "Meshes" include("meshes/00-runtests.jl")
    @testset "Mathematics" include("mathematics/00-runtests.jl")    
    @testset "Utilities" include("utilities/00-runtests.jl")    
    @testset "Plots" include("plots/00-runtests.jl")    
    @testset "Gmsh" include("gmsh/00-runtests.jl")    
end

end
