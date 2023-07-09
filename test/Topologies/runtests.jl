module TopologiesTests

using Test

@testset "ConnectivityList" begin include("ConnectivityListTests.jl") end

@testset "Topology" begin include("TopologyTests.jl") end

end # module

