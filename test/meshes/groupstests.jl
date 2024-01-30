using Test

using MMJMesh
using MMJMesh.Meshes


# Test struct
struct Foo{T}
    index::Int
end
MMJMesh.Meshes.index(f::Foo) = f.index


# -------------------------------------------------------------------------------------------------
# Group
# -------------------------------------------------------------------------------------------------

# Set up
f11 = Foo{1}(1)
f21 = Foo{2}(1)
f13 = Foo{1}(3)
f31 = Foo{3}(1)
f33 = Foo{3}(3)

# Group of general type Foo
g = Group{Foo}(1:2);
@test 99 ∉ g
@test f11 ∈ g
@test f21 ∈ g
@test f13 ∉ g

# Group of specialized type Foo{1}
g = Group{Foo{1}}(1:2);
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
    gc.entries[:zab] = Group{Foo{1}}(1:2)
    return nothing
end

gc[:foo] = Group{Foo{1}}(1:2)
gc[:bar] = Group{Foo{2}}(1:2)
addrecipe!(gc, :baz, () -> Group{Foo{3}}(1:length(a)))  
addrecipe!(gc, :zab, tr)  

@test length(gc[:foo]) == 2
@test f11 ∈ gc[:foo]
@test f21 ∉ gc[:foo]
@test f13 ∉ gc[:foo]

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
