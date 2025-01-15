using Test
using IntervalSets
using DomainSets: ×

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Utilities
using MMJMesh.Geometries

using MMJMesh.Meshes: geometrytype


# Set up
coords = [0.0 1.0 2.0 0.1 0.9 1.9; 0.0 0.1 0.0 0.9 1.0 0.9]
elts = [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
m = Mesh(coords, elts, 2)


# -------------------------------------------------------------------------------------------------
# Dimensions of types
# -------------------------------------------------------------------------------------------------

@test pdim(Node) == 0
@test pdim(Edge) == 1
@test pdim(Face) == 2
@test pdim(Solid) == 3

@test gdim(Node{1}) == 1
@test gdim(Edge{2}) == 2
@test gdim(Face{3}) == 3
@test gdim(Solid{3}) == 3

@test nnodes(Node) == 0
@test nnodes(Edge{2,3}) == 3
@test nnodes(Face{3,4}) == 4


# -------------------------------------------------------------------------------------------------
# Node
# -------------------------------------------------------------------------------------------------

n5 = node(m, 5)
@test pdim(n5) == 0
@test gdim(n5) == 2
@test nedges(n5) == 3
@test coordinates(n5) == [0.9, 1.0]
@test coordinate(n5, 1) == 0.9
@test coordinate(n5, 2) == 1.0


# -------------------------------------------------------------------------------------------------
# Edge
# -------------------------------------------------------------------------------------------------

e2 = edge(m, 2)
@test pdim(e2) == 1
@test gdim(e2) == 2
@test nnodes(e2) == 2
@test length(e2) == hypot(0.1, 0.9)

g = geometry(e2)
p = parametrization(g)
@test p(-1) == coordinates(e2, 1)
@test p(1) == coordinates(e2, 2)


# -------------------------------------------------------------------------------------------------
# Face
# -------------------------------------------------------------------------------------------------

f1 = face(m, 1)
@test pdim(f1) == 2
@test gdim(f1) == 2
@test nnodes(f1) == 4
@test nodeindices(f1) == [1, 2, 5, 4]
@test edgeindices(f1) == [1, 2, 3, 4]
@test faceindices(f1) == [3]
@test coordinates(f1) == [0.0 1.0 0.9 0.1; 0.0 0.1 1.0 0.9]
@test coordinates(f1, 1) == [0.0; 0.0]
@test coordinates(f1, 2) == [1.0; 0.1]
@test coordinates(f1, 4) == [0.1; 0.9]

f3 = face(m, 3)
@test nnodes(f3) == 3
@test nodeindices(f3) == [2, 6, 5]
@test coordinates(f3) == [1.0 1.9 0.9; 0.1 0.9 1.0]
@test coordinates(f3, 1) == [1.0, 0.1]
@test coordinates(f3, 2) == [1.9, 0.9]
@test coordinates(f3, 3) == [0.9, 1.0]


# -------------------------------------------------------------------------------------------------
# Geometry
# -------------------------------------------------------------------------------------------------

m = Mesh((0 .. 4) × (0 .. 2), 4, 1)
n = node(m, 1)
f = face(m, 3)
phi = parametrization(geometry(f))

@test geometry(n) isa Point
@test geometrytype(f) == Box
@test geometry(f) isa Box

@test phi(-1, -1) == coordinates(f, 1)
@test phi(1, -1) == coordinates(f, 2)
@test phi(1, 1) == coordinates(f, 3)
@test phi(-1, 1) == coordinates(f, 4)


# -------------------------------------------------------------------------------------------------
# Access to parts
# -------------------------------------------------------------------------------------------------

m = Mesh((0 .. 1) × (0 .. 1), 3, 3)

e = element(m, 1)
@test node(e, 1) === node(m, 1)
@test node(e, 2) === node(m, 2)
@test node(e, 3) === node(m, 6)
@test node(e, 4) === node(m, 5)
@test edge(e, 1) == edge(m, 1)

e = element(m, 9)
@test node(e, 1) === node(m, 11)
@test edge(e, 3) === edge(m, 24)


# -------------------------------------------------------------------------------------------------
# Access multiple coordinates
# -------------------------------------------------------------------------------------------------

m = Mesh(meshpath("complex-g1.msh"))
@test coordinates(m, [3, 9]) == stack([coordinates(m, 3), coordinates(m, 9)])
@test coordinates(face(m, 1), [1, 3]) == coordinates(m, [329, 176])
