module Geometries

# Exports
## geometry.jl
export Geometry
## point.jl
export Point
export coordinates

# Modules needed by this module
using StaticArrays
using MMJMesh.MMJBase

# Parts
include("detail.jl")
include("point.jl")
include("geometry.jl")

end