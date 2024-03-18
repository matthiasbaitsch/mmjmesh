"""
    GeometricObject{DT, DG}

A geometric object of parametric dimension `DT` embedded in `DG` dimensional space.
Here, a geometric object is understood as a connected set of points.
"""
abstract type GeometricObject{DT,DG} end

pdim(o::GeometricObject) = pdim(typeof(o))
gdim(o::GeometricObject) = gdim(typeof(o))
gdim(::Type{<:GeometricObject}) = @notimplemented
pdim(::Type{<:GeometricObject}) = @notimplemented
center(::GeometricObject) = @notimplemented
measure(::GeometricObject) = @notimplemented
boundingbox(::GeometricObject) = @notimplemented
Base.in(::GeometricObject, x; eps::Float64=0.0) = @notimplemented


"""
    GeometricObjectP{DT,DG}

A `DT,DG` dimensional geometric object with a parametrization. For such an object, the geometry is the 
image of a simple `DT`-dimensional reference domain into `DG` dimensions.
"""
abstract type GeometricObjectP{DT,DG} <: GeometricObject{DT,DG} end

parametrization(::GeometricObjectP) = @notimplemented


"""
    GeometricObjectI{NP,DT,DG}

A `DT,DG` dimensional geometric object with a parametrization which interpolates `NP` points. An
important class of such objects are the geometries of typical isoparametric finite elements.
"""
struct GeometricObjectI{DT,DG,NP} <: GeometricObjectP{DT,DG}

end