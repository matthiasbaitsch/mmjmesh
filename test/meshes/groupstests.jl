using Test

using MMJMesh
using MMJMesh.Meshes


# Test struct
struct Bar{T}
    index::Int
end
MMJMesh.Meshes.index(f::Bar) = f.index


# -------------------------------------------------------------------------------------------------
# Group
# -------------------------------------------------------------------------------------------------

# Set up
f11 = Bar{1}(1)
f21 = Bar{2}(1)
f13 = Bar{1}(3)
f31 = Bar{3}(1)
f33 = Bar{3}(3)

# Group of general type Bar
g = Group{Bar}(1:2);
@test 99 ∉ g
@test f11 ∈ g
@test f21 ∈ g
@test f13 ∉ g

# Group of specialized type Bar{1}
g = Group{Bar{1}}(1:2);
@test 99 ∉ g
@test f11 ∈ g
@test f21 ∉ g
@test f13 ∉ g


# -------------------------------------------------------------------------------------------------
# GroupCollection
# -------------------------------------------------------------------------------------------------
a = [1, 2, 6]
gc = GroupCollection()

function tr()
    global gc
    gc.entries[:zab] = Group{Bar{1}}(1:2)
    return nothing
end

gc[:Bar] = Group{Bar{1}}(1:2)
gc[:bar] = Group{Bar{2}}(1:2)
addrecipe!(gc, :baz, () -> Group{Bar{3}}(1:length(a)))  
addrecipe!(gc, :zab, tr)  

@test length(gc[:Bar]) == 2
@test f11 ∈ gc[:Bar]
@test f21 ∉ gc[:Bar]
@test f13 ∉ gc[:Bar]

@test length(gc[:bar]) == 2
@test f11 ∉ gc[:bar]
@test f21 ∈ gc[:bar]
@test f13 ∉ gc[:bar]

@test isnothing(gc.entries[:baz])
@test length(gc[:baz]) == 3
@test f11 ∉ gc[:baz]
@test f31 ∈ gc[:baz]
@test !isnothing(gc.entries[:baz])
clearcache!(gc)
@test isnothing(gc.entries[:baz])
push!(a, 99)
@test length(gc[:baz]) == 4

@test length(gc[:zab]) == 2
@test f11 ∈ gc[:zab]
