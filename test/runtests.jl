module MMJMeshRunTests

using Test

@testset "Topologies" include("topologies/runtests.jl")

@testset "Geometries" include("geometries/runtests.jl")

end
