using Test
using MMJMesh.Meshes

# -------------------------------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------------------------------

# Test mesh
m = Mesh(:quadtri)
definegroup!(:bar, m, 2, [2])

# Assign data
setdata!(m, :foo, 1)
setdata!(group(m, :bar), :foo, 2)
setdata!(group(m, :elements), :foo, 3)
setdata!(node(m, 1), :foo, 4)
setdata!(element(m, 1), :foo, 5)

# Test data(...)
@test data(m, :foo) == 1
@test data(group(m, :bar), :foo) == 2
@test data(group(m, :elements), :foo) == 3
@test data(node(m, 1), :foo) == 4
@test data(element(m, 1), :foo) == 5
@test data(element(m, 2), :foo) == 2
@test data(element(m, 3), :foo) == 3
@test data(edge(m, 1), :foo) == 1
@test data(edge(m, 1), :baz) === nothing

# Test hasdata(...)
@test hasdata(edge(m, 1), :foo)
@test !hasdata(edge(m, 1), :baz)

# Test reassign data
setdata!(m, :foo, 2)
setdata!(group(m, :bar), :foo, 3)
setdata!(group(m, :elements), :foo, 4)
setdata!(node(m, 1), :foo, 5)
setdata!(element(m, 1), :foo, 6)
@test data(m, :foo) == 2
@test data(group(m, :bar), :foo) == 3
@test data(group(m, :elements), :foo) == 4
@test data(node(m, 1), :foo) == 5
@test data(element(m, 1), :foo) == 6
@test data(element(m, 2), :foo) == 3
@test data(element(m, 3), :foo) == 4
@test data(edge(m, 1), :foo) == 2
@test data(edge(m, 1), :baz) === nothing

## XXX

# -------------------------------------------------------------------------------------------------
# Test array data
# -------------------------------------------------------------------------------------------------

# Test mesh
m = Mesh(:quadtri)

# Assign data
setdata!(m, :foo, 99)
setdata!(m, :bar, [4, 1, 8])

# Test data(...)
@test data(element(m, 1), :foo) == 99
@test data(element(m, 1), :bar) == 4
@test data(element(m, 2), :bar) == 1
@test data(element(m, 3), :bar) == 8


