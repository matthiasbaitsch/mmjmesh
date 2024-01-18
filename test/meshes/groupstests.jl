using Test

using MMJMesh
using MMJMesh.Meshes


# Set up
struct Foo{T}
    index::Int
end
f11 = Foo{1}(1)
f21 = Foo{2}(1)
f13 = Foo{1}(3)


# Group of general type Foo
g = EntityGroup{Foo}(1:2);
@test 99 ∉ g
@test f11 ∈ g
@test f21 ∈ g
@test f13 ∉ g


# Group of specialized type Foo{1}
g = EntityGroup{Foo{1}}(1:2);
@test 99 ∉ g
@test f11 ∈ g
@test f21 ∉ g
@test f13 ∉ g
