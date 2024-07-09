module Utilities

# Packages needed
using DomainSets
using IntervalSets

using DomainSets: ×

# Modules needed by this module
using MMJMesh.Meshes
using MMJMesh.Topologies
using MMJMesh.Geometries

# Exports
## generatemeshes.jl
export Meshtype, QUADRANGLE, TRIANGLE
export makemeshoninterval, makemeshonrectangle

# Parts
include("generatemeshes.jl")

end
