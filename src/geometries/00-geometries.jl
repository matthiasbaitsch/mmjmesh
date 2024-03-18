module Geometries

# Modules needed by this module
using StaticArrays
using LinearAlgebra

using MMJMesh
using MMJMesh.MMJBase
using MMJMesh.Mathematics

# Functions extended by this module
import Base.in
import Base.isequal

import MMJMesh.MMJBase: gdim, pdim

# Parts
include("detail.jl")
include("geometricobject.jl")
include("point.jl")
include("box.jl")
include("geometry.jl")

# Exports
## geometricobject.jl
export GeometricObject, GeometricObjectP, GeometricObjectI,
       parametrization,
       center, measure, boundingbox
## geometry.jl
export Geometry
## point.jl
export Point
## box.jl
export Box, diagonal, sides
## various
# export coordinates

end

