module MeshesTests

using Test

@testset "Common" include("datatests.jl")
@testset "Common" include("groupstests.jl")
@testset "Meshes" include("meshtests.jl")
@testset "Meshes" include("meshentitiestests.jl")
@testset "Common" include("entitydatatests.jl")
@testset "Common" include("entitygroupstests.jl")

end

