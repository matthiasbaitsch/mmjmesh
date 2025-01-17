module GeometriesTests

using Test

@testset "Point" include("pointtests.jl")
@testset "Box" include("boxtests.jl")
@testset "GeometricObject" include("geometricobjecttests.jl")
@testset "Geometry" include("geometrytests.jl")
@testset "Lines" include("linestests.jl")

end

