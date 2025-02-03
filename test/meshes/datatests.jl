using Test
using MMJMesh.Meshes

# Create test mesh
coords = [0.0 1.0 2.0 0.1 0.9 1.9; 0.0 0.1 0.0 0.9 1.0 0.9]
elts = [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
m = Mesh(coords, elts, 2)
definegroup!(m, 2, :bar, [2])

# Assign data
setdata!(m, :foo, 1)
setdata!(group(m, :bar), :foo, 2)
setdata!(group(m, :elements), :foo, 3)
setdata!(node(m, 1), :foo, 4)
setdata!(element(m, 1), :foo, 5)

# Test direct access
@test data(m, :foo) == 1
@test data(group(m, :bar), :foo) == 2
@test data(group(m, :elements), :foo) == 3
@test data(node(m, 1), :foo) == 4
@test data(element(m, 1), :foo) == 5
@test data(element(m, 2), :foo) == 2
@test data(element(m, 3), :foo) == 3
@test data(edge(m, 1), :foo) == 1
@test data(edge(m, 1), :baz) === nothing

@test hasdata(edge(m, 1), :foo)
@test !hasdata(edge(m, 1), :baz)

# Test array access
coords = [0.0 1.0 2.0 0.1 0.9 1.9; 0.0 0.1 0.0 0.9 1.0 0.9]
elts = [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
m = Mesh(coords, elts, 2)
setdata!(m, :foo, 99)
setdata!(m, :bar, [4, 1, 8])

@test data(element(m, 1), :foo) == 99
@test data(element(m, 1), :bar) == 4
@test data(element(m, 2), :bar) == 1
@test data(element(m, 3), :bar) == 8
