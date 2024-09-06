using Test
using MMJMesh
using MMJMesh.Meshes

import MMJMesh.Meshes: MeshData

# -------------------------------------------------------------------------------------------------
# Global data
# -------------------------------------------------------------------------------------------------
md = MeshData{Int}()
md[:foo] = 42

@test md[:foo] == 42
@test md[:foo, 22] == 42


# -------------------------------------------------------------------------------------------------
# Actual mapping using a callable object
# -------------------------------------------------------------------------------------------------

mutable struct Foo
    a::Int
    b::Vector{Int}
    value::Int
end

(f::Foo)(p1::Int, p2::Int) = p1 == f.a && p2 in f.b ? f.value : nothing

md.mappings[:bar] = Foo(2, [1, 2, 33], 43)

@test md[:bar, 2, 1] == 43
@test md[:bar, 2, 2] == 43
@test md[:bar, 2, 33] == 43
@test isnothing(md[:bar, 3, 33])

md.mappings[:bar].value = 44
@test md[:bar, 2, 1] == 44

