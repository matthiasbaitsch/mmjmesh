using Test

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Utilities


# -------------------------------------------------------------------------------------------------
# Test for 1D mesh
# -------------------------------------------------------------------------------------------------
m = makemeshoninterval(0.0, 8.0, 7)
e5 = element(m, 5)

e5.data[:baz] = 66
@test e5.data[:baz] == 66

m.data[:foo] = 88
@test e5.data[:foo] == 88

m.data[:bar, :edges] = 77
@test e5.data[:bar] == 77


# -------------------------------------------------------------------------------------------------
# Test for 2D mesh
# -------------------------------------------------------------------------------------------------

# Prepare test data
a = 5
m = makemeshonrectangle(9.0, 4.5, 2a, a)
definegroup!(m, 2, :g1, [1, 2, 3, 13, 14, 7, 15])
definegroup!(m, 2, :g2, [1, 45, 56])

e1 = edge(m, 1)
f1 = element(m, 1)
f13 = element(m, 13)
f42 = element(m, 42)
f45 = element(m, 45)

# Associate with everything
m.data[:d1] = 42
@test m.data[:d1] == 42
@test f42.data[:d1] == 42

# Associate with group
m.data[:d2, :g1] = 13
m.data[:d2, :g2] = 81
@test f1.data[:d2] == 81
@test f13.data[:d2] == 13
@test isnothing(f42.data[:d2])
@test f45.data[:d2] == 81

# Associate with all faces
m.data[:d3, :g1] = 91
m.data[:d3, :faces] = 92
@test f13.data[:d3] == 91
@test f42.data[:d3] == 92
@test isnothing(e1.data[:d3])
@test_throws InexactError m.data[:d3, :g2] = 91.1

# Different functions with one name
m.data[:d4, :g1] = sin
m.data[:d4, :g2] = cos
@test f13.data[:d4](π / 2) ≈ 1.0
@test f45.data[:d4](π) ≈ -1.0

# Associate with mesh entities
e1.data[:d5] = 11
f13.data[:d5] = 12
@test e1.data[:d5] == 11
@test f13.data[:d5] == 12
@test isnothing(f42.data[:d5])

# Associate with all nodes
m.data[:d6, 0] = collect(1:nnodes(m))

for i ∈ 1:nnodes(m)
    @test node(m, i).data[:d6] == i
end

# Different functions with one name
e1.data[:d7] = sin
f13.data[:d7] = cos
@test e1.data[:d7](π / 2) ≈ 1.0
@test f13.data[:d7](π) ≈ -1.0
