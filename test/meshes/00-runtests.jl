module MeshesTests

using Test

@testset "Data" include("datatests.jl")
@testset "Group" include("groupstests.jl")
@testset "Mesh" include("meshtests.jl")
@testset "MeshEntity" include("meshentitiestests.jl")
@testset "EntityData" include("entitydatatests.jl")
@testset "MeshEntityGroup" include("entitygroupstests.jl")

end

