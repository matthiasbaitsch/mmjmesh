module TopologiesTests

using Test

@testset "ConnectivityList" include("connectivitylisttests.jl")

@testset "Topology" include("topologytests.jl")

end

