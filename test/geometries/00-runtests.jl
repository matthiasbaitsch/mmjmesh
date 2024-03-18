module GeometriesTests

using Test

@testset "Point" include("pointtests.jl")
@testset "Box" include("boxtests.jl")
@testset "Geometry" include("geometrytests.jl")

end

