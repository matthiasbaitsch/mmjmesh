"""
    Box(p1, p2)

Box in d-dimensional space defined by two points `p1` and `p2`.
"""
struct Box{D} <: GeometricObject{D,D}
    min::Point{D}
    max::Point{D}

    function Box{D}(min, max) where {D}
        @assert min ≤ max "`min` must be less than or equal to `max`"
        new(min, max)
    end
end

Box(p1::Point{D}, p2::Point{D}) where {D} = Box{D}(p1, p2)
Box(p1::Tuple, p2::Tuple) = Box(Point(p1...), Point(p2...))

MMJMesh.gdim(::Type{<:Box{D}}) where {D} = D
MMJMesh.pdim(::Type{<:Box{D}}) where {D} = D

center(b::Box) = Point((coordinates(b.max) + coordinates(b.min)) / 2)
diagonal(b::Box) = norm(b.max - b.min)
sides(b::Box) = Tuple(b.max - b.min)

Base.minimum(b::Box) = b.min
Base.maximum(b::Box) = b.max
Base.extrema(b::Box) = b.min, b.max
Base.isapprox(b1::Box, b2::Box) = b1.min ≈ b2.min && b1.max ≈ b2.max
Base.in(p::Point, b::Box) = b.min ≤ p ≤ b.max
