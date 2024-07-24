"""
    GeometricObject{DT, DG}

A geometric object of parametric dimension `DT` embedded in `DG` dimensional space.
Here, a geometric object is understood as a connected set of points.
"""
abstract type GeometricObject{DT,DG} end

MMJMesh.pdim(o::GeometricObject) = pdim(typeof(o))
MMJMesh.gdim(o::GeometricObject) = gdim(typeof(o))
MMJMesh.gdim(::Type{<:GeometricObject}) = @notimplemented
MMJMesh.pdim(::Type{<:GeometricObject}) = @notimplemented

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
    GeometricObjectI{DT,DG,NP}

A `DT,DG` dimensional geometric object with a parametrization which interpolates `NP` points. An
important class of such objects are the geometries of isoparametric finite elements.
"""
struct GeometricObjectI{DT,DG,NP} <: GeometricObjectP{DT,DG}
    points::SMatrix{DG,NP,Float64}

    GeometricObjectI{DT}(points::AbstractMatrix) where {DT} =
        new{DT,size(points, 1),size(points, 2)}(points)
end

parametrization(o::GeometricObjectI{DT}) where {DT} =
    Interpolation(bases[(DT, size(o.points, 2))], o.points)

"""
Lookup `(DT, NN) → basis` where `DT` is the topological dimension and `NN` the number of nodes
of the geometric entity.
"""
const bases = Dict(
    (1, 2) => MappingFromComponents(nodalbasis(makeelement(:lagrange, IHat, k=1))...),
    (2, 4) => MappingFromComponents(nodalbasis(makeelement(:lagrange, QHat, k=1))...)
)

