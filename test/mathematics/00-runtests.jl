module MathematicsTests

using Test

include("Validate.jl")

@testset "Domains"          include("domainstests.jl")
@testset "Mappings"         include("mappingstests.jl")
@testset "MPolynomials"     include("mpolynomials2tests.jl")
@testset "Piecewise"        include("piecewisetests.jl")
@testset "Curves"           include("curvestests.jl")
@testset "Interpolations"   include("interpolationstests.jl")
@testset "Finiteelements"   include("finiteelementstests.jl")

end

