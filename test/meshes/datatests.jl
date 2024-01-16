using Test
using MMJMesh
using MMJMesh.Meshes

import MMJMesh.Meshes: Data, addmapping!, setbase!


# -------------------------------------------------------------------------------------------------
# Global data
# -------------------------------------------------------------------------------------------------
md = Data(Int)
md[:foo] = 42
@test md[:foo] == 42


# -------------------------------------------------------------------------------------------------
# Mapping based on dimension and set
# -------------------------------------------------------------------------------------------------
function makesetmapping(dim::Int, set::Set, value::Any)
    function mapping(base::Int, d::Int, i::Int)
        if d == dim && i in set
            return base + value
        end
        return nothing
    end
end

setbase!(md, 11)
addmapping!(md, :bar, makesetmapping(2, Set([1, 2, 33]), 43))

@test md[:bar, 2, 1] == 11 + 43
@test md[:bar, 2, 2] == 11 + 43
@test md[:bar, 2, 33] == 11 + 43
@test isnothing(md[:bar, 2, 3])
@test isnothing(md[:bar, 1, 1])
