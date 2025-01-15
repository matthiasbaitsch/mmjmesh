module MathematicsTests

using Test

using MMJMesh.Tests
using MMJMesh.Mathematics

@testset "Mappings" begin
    include("domainstests.jl")
    include("mappingstests.jl")
    include("mpolynomialstests.jl")
    include("piecewisetests.jl")
    include("curvestests.jl")
    include("interpolationstests.jl")
    include("finiteelementstests.jl")
end

end

