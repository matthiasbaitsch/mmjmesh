module UtilitiesTests

using Test

@testset "generatemeshes" include("generatemeshestests.jl")
@testset "SeqIntSet" include("seqintsettests.jl")

end

