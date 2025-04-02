module Utilities


# Packages needed
import DomainSets
import IntervalSets

using DomainSets: ×
using IntervalSets: (..)


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
