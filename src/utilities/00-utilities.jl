module Utilities


# Packages needed
import DomainSets
import IntervalSets

# Modules needed by this module
using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Topologies
using MMJMesh.Geometries
using MMJMesh.Mathematics


# Exports

## generatemeshes.jl
export Meshtype, QUADRANGLE, TRIANGLE, CRISSCROSS

# Parts
include("generatemeshes.jl")

end
