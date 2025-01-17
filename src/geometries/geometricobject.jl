# -------------------------------------------------------------------------------------------------
# Basic geometric object
# -------------------------------------------------------------------------------------------------

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
Base.in(::RealVec, ::GeometricObject; atol::Real=0.0) = @notimplemented


# -------------------------------------------------------------------------------------------------
# Geometric object with parametrization
# -------------------------------------------------------------------------------------------------

"""
    GeometricObjectP{DT,DG}

A `DT,DG` dimensional geometric object with a parametrization. For such an object, the geometry is the 
image of a simple `DT`-dimensional reference domain into `DG` dimensions.
"""
abstract type GeometricObjectP{DT,DG} <: GeometricObject{DT,DG} end

"""
    parametrization(o)

Returns the map `F` such that the geometric object `o` is the image of `F`.
"""
parametrization(::GeometricObjectP) = @notimplemented

"""
    parameterof(o, p, atol=1e-12)

Returns the parameter of point `p` on the one-dimensional geometric object o. The functions returns
`NaN` if `p` is further than `atol` away from `o`.
"""
parameterof(o::GeometricObjectP{1}, p; atol::Real=1e-12) = @notimplemented

Base.in(p::RealVec, o::GeometricObject{1}; atol::Real=0.0) = !isnan(parameterof(o, p, atol=atol))

# -------------------------------------------------------------------------------------------------
# Geometric object from point interpolation
# -------------------------------------------------------------------------------------------------

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

measure(e::GeometricObjectI{1,DG,2}) where {DG} = norm(diff(e.points, dims=2))
measure(f::GeometricObjectI{2,2,NP}) where {NP} =
    0.5 * sum([det(f.points[:, [i, i % NP + 1]]) for i = 1:NP])

"""
Lookup `(DT, NN) → basis` where `DT` is the topological dimension and `NN` the number of nodes
of the geometric entity.
"""
const bases = Dict(
    (1, 2) => MappingFromComponents(nodalbasis(makeelement(:lagrange, IHat, k=1))...),
    (2, 4) => MappingFromComponents(nodalbasis(makeelement(:lagrange, QHat, k=1))...)
)

