module Utilities

# Packages needed
using DomainSets
using IntervalSets

using DomainSets: ×

# Modules needed by this module
using MMJMesh.Meshes
using MMJMesh.Topologies
using MMJMesh.Geometries
using MMJMesh.Mathematics

# Exports
## generatemeshes.jl
export Meshtype, QUADRANGLE, TRIANGLE

# Parts
include("generatemeshes.jl")

end
