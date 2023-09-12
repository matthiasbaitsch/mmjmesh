module Geometries


export Geometry
export coordinates

export Point

using StaticArrays
using MMJMesh.MMJBase


include("details.jl")

include("point.jl")

include("geometry.jl")

end