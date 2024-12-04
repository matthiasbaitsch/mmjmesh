module MeshesTests

using Test

@testset "Data" include("datatests.jl")
@testset "Group" include("groupstests.jl")
@testset "Mesh" include("meshtests.jl")
@testset "MeshEntity" include("meshentitiestests.jl")
@testset "MeshEntityGroup" include("entitygroupstests.jl")
@testset "Build" include("buildtests.jl")

end

