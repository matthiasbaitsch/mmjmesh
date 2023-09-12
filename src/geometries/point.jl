"""
    Point(x1, ..., xD)

    A point in `D`-dimensional space. Based on `point.jl` from Meshes.jl.
"""
struct Point{D}
    coordinates::SVector{D,Float64}
end

Point{D}(coords...) where D = Point{D}(coords)

Point(a1::Real, a2::Real) = Point{2}(a1, a2)
Point(a1::Real, a2::Real, a3::Real) = Point{3}(a1, a2, a3)

pdim(A::Point) = 0
gdim(::Point{D}) where D = D
center(A::Point) = A
measure(::Point) = 0.0
boundingbox(::Point) = nothing
coordinates(A::Point) = A.coordinates

Base.:(-)(A::Point, B::Point) = A.coordinates - B.coordinates
Base.:(==)(A::Point, B::Point) = A.coordinates == B.coordinates
Base.in(A::Point, B::Point) = A == B
Base.isapprox(A::Point{D}, B::Point{D}; atol=atol(Float64), kwargs...) where D = isapprox(A.coords, B.coords; atol, kwargs...)

function Base.show(io::IO, point::Point)
    print(io, "Point$(Tuple(point.coordinates))")
end
