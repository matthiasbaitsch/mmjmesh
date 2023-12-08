"""
    Point(x1, ..., xD)

    A point in `D`-dimensional space. Based on `point.jl` from Meshes.jl.
"""
struct Point{D} <: GeometricObject{0, D}
    coordinates::SVector{D,Float64}
end

Point{D}(coords...) where {D} = Point{D}(coords)
Point(a1::Number, a2::Number) = Point{2}(a1, a2)
Point(a1::Number, a2::Number, a3::Number) = Point{3}(a1, a2, a3)

pdim(A::Point) = 0
gdim(::Point{D}) where {D} = D
center(A::Point) = A
measure(::Point) = 0.0
boundingbox(::Point) = nothing
coordinates(A::Point) = A.coordinates

Base.:(-)(A::Point, B::Point) = A.coordinates - B.coordinates
Base.:(==)(A::Point, B::Point) = A.coordinates == B.coordinates
Base.in(A::Point, B::Point) = A == B
Base.isapprox(A::Point, B::Point; atol=atol(Float64), kwargs...) = isapprox(A.coordinates, B.coordinates; atol, kwargs...)

function Base.show(io::IO, point::Point)
    print(io, "Point$(Tuple(point.coordinates))")
end
