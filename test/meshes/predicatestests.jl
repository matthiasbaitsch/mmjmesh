using Test

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Geometries

# Elementary example
coords = [0 1 2 0 1 2; 0 0 0 1 1 1]
elts = [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
m = Mesh(coords, elts, 2)

p = on(HLine(1))
@test isnan(p(node(m, 1)))
@test p(node(m, 4)) == 0
@test p(node(m, 5)) == 1
@test isnan(p(edge(m, 2)))
@test p(edge(m, 3)) == 0.5


# Complex example
meshpath(m) = joinpath(@__DIR__(), "../../data/gmsh", m)

m = Mesh(meshpath("complex-g1.msh"))

# Nodes
@test nodeindex(m, on(Point(0, 0))) == 5
@test nodeindices(m, on(Point(0, 0))) == [5]
@test nodeindices(m, on(Segment([0, 0], [0.5, 0]))) == [5, 46, 47, 48, 49, 50, 51, 52, 53, 54, 6]
@test nodeindices(m, on(HLine(0))) ==
      [5, 46, 47, 48, 49, 50, 51, 52, 53, 54, 6, 9, 89, 90, 91, 92, 93, 94, 95, 96, 97, 10]
@test nodeindices(m, on(VLine(0))) ==
      [5, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 8]

# Edges
@test edgeindices(m, on(VLine(0))) == 
      [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21]



 