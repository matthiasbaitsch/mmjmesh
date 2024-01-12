module MeshesTests

using Test

@testset "Meshes" include("meshestests.jl")
@testset "Common" include("entitygroupstests.jl")

end

