module MMJBaseTests

using Test

@testset "SeqIntSet" include("seqintsettests.jl")
@testset "Functions" include("functionstests.jl")

end

