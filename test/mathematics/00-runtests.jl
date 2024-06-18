module MathematicsTests

using Test

@testset "Mappings" begin
    include("mappingstests.jl")
    include("mpolynomialstests.jl")
end

end

