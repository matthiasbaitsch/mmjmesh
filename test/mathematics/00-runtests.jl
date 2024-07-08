module MathematicsTests

using Test

using MMJMesh.Mathematics
include("validatemappings.jl")

@testset "Mappings" begin
    include("mappingstests.jl")
    include("mpolynomialstests.jl")
    include("fefunctionstests.jl")
end

end

