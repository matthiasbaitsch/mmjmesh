module Geometries

# Modules needed by this module
import DomainSets
using StaticArrays
using LinearAlgebra

using MMJMesh
using MMJMesh.MMJBase
using MMJMesh.Mathematics
import MMJMesh: coordinate, coordinates

# Exports
## geometricobject.jl
export GeometricObject, GeometricObjectP, GeometricObjectI,
    parametrization, parameterof,
    center, measure, boundingbox

## geometry.jl
export Geometry

## point.jl
export Point

## box.jl
export Box, diagonal, sides

## lines.jl
export HLine, VLine, Segment


# Parts
include("detail.jl")
include("geometricobject.jl")
include("point.jl")
include("box.jl")
include("geometry.jl")
include("lines.jl")

end

