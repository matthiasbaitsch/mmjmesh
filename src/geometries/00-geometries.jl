module Geometries

# Modules needed by this module
using StaticArrays
using MMJMesh
using MMJMesh.MMJBase

# Functions extended by this module
import Base.in
import Base.isequal

# Parts
include("detail.jl")
include("geometricobject.jl")
include("point.jl")
include("geometry.jl")

# Exports
## geometricobject.jl
export GeometricObject
## geometry.jl
export Geometry
## point.jl
export Point

end