module MeshesTests

using Test

@testset "Data" include("datatests.jl")
@testset "Group" include("groupstests.jl")
@testset "Mesh" include("meshestests.jl")
@testset "MeshEntity" include("meshentitiestests.jl")
@testset "MeshEntityGroup" include("meshentitygroupstests.jl")
@testset "Build" include("buildtests.jl")
@testset "Predicates" include("predicatestests.jl")

end

