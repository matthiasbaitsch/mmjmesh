using Test

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Geometries


meshpath(m) = joinpath(@__DIR__(), "../../data/gmsh", m)

m = Mesh(meshpath("complex-g1.msh"))

@test nodeindex(m, on(Point(0, 0))) == 5
@test nodeindices(m, on(Point(0, 0))) == [5]
@test nodeindices(m, on(Segment([0, 0], [0.5, 0]))) == [5, 46, 47, 48, 49, 50, 51, 52, 53, 54, 6]
@test nodeindices(m, on(HLine(0))) ==
      [5, 46, 47, 48, 49, 50, 51, 52, 53, 54, 6, 9, 89, 90, 91, 92, 93, 94, 95, 96, 97, 10]
@test nodeindices(m, on(VLine(0))) ==
      [5, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 8]