using Test
using MMJMesh
using MMJMesh.Meshes


import MMJMesh.Meshes: MeshData, addmapping!

# Global data
md = MeshData(Int)
md[:foo] = 42

@test md[:foo] == 42


# Mapping based on dimension and set
function makesetmapping(dim::Int, set::Set, value::Any)
    function mapping(p...)
        (_, d, i) = p
        if d == dim && i in set
            return value
        end
        return nothing
    end
end
addmapping!(md, :bar, makesetmapping(2, Set([1, 2, 33]), 43))

@test md[:bar, 2, 1] == 43
@test md[:bar, 2, 2] == 43
@test md[:bar, 2, 33] == 43
@test isnothing(md[:bar, 2, 3])
@test isnothing(md[:bar, 1, 1])
