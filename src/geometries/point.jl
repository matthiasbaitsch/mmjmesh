"""
    Point(x1, ..., xD)

    p1 point in `D`-dimensional space. Based on `point.jl` from Meshes.jl.
"""
struct Point{D} <: GeometricObject{0,D}
    coordinates::SVector{D,Float64}
end

Point{D}(coords...) where {D} = Point{D}(coords)
Point(a1::Number, a2::Number) = Point{2}(a1, a2)
Point(a1::Number, a2::Number, a3::Number) = Point{3}(a1, a2, a3)

MMJMesh.pdim(::Type{<:Point}) = 0
MMJMesh.gdim(::Type{<:Point{D}}) where {D} = D

coordinates(p::Point) = p.coordinates
center(p::Point) = p
measure(::Point) = 0.0
boundingbox(::Point) = nothing

Base.:(-)(p1::Point, p2::Point) = p1.coordinates - p2.coordinates
Base.:(==)(p1::Point, p2::Point) = p1.coordinates == p2.coordinates
Base.in(p1::Point, p2::Point) = p1 == p2
Base.isapprox(p1::Point, p2::Point; atol=atol(Float64), kwargs...) = isapprox(p1.coordinates, p2.coordinates; atol, kwargs...)

Base.:(≤)(p1::Point{D}, p2::Point{D}) where {D} = all(≥(zero(D)), p2 - p1)
Base.:(≥)(p1::Point{D}, p2::Point{D}) where {D} = all(≥(zero(D)), p1 - p2)
Base.:(<)(p1::Point{D}, p2::Point{D}) where {D} = all(>(zero(D)), p2 - p1)
Base.:(>)(p1::Point{D}, p2::Point{D}) where {D} = all(>(zero(D)), p1 - p2)

function Base.show(io::IO, point::Point)
    print(io, "Point$(Tuple(point.coordinates))")
end
