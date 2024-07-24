module MathematicsTests

using Test

using MMJMesh.Mathematics
include("validatemappings.jl")

@testset "Mappings" begin
    include("domainstests.jl")
    include("mappingstests.jl")
    include("mpolynomialstests.jl")
    include("fefunctionstests.jl")
    include("finiteelementstests.jl")
end

end

