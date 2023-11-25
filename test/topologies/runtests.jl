module TopologiesTests

using Test

@testset "ConnectivityList" include("connectivitylisttests.jl")

@testset "Topology2D" include("topologytests2d.jl")

end

