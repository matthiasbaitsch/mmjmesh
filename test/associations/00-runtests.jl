module AssociationsTests

using Test

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Utilities
using MMJMesh.Associations

m = makemeshoninterval(0.0, 1.2, 10)

m.data[:test] = 1
@test m.data[:test] == 1

end

