"""
    Box(p1, p2)

Box in d-dimensional space defined by two points `p1` and `p2`.
"""
struct Box{D} <: GeometricObjectP{D,D}
    min::Point{D}
    max::Point{D}

    function Box{D}(min, max) where {D}
        @assert min ≤ max "`min` must be less than or equal to `max`"
        new(min, max)
    end
end

Box(p1::Point{D}, p2::Point{D}) where {D} = Box{D}(p1, p2)
Box(p1::AbstractArray, p2::AbstractArray) = Box(Point(p1...), Point(p2...))
Box(r::DomainSets.Rectangle) = Box(DomainSets.leftendpoint(r), DomainSets.rightendpoint(r))

MMJMesh.gdim(::Type{<:Box{D}}) where {D} = D
MMJMesh.pdim(::Type{<:Box{D}}) where {D} = D

center(b::Box) = Point((coordinates(b.max) + coordinates(b.min)) / 2)
diagonal(b::Box) = norm(b.max - b.min)
sides(b::Box) = b.max - b.min

parametrization(b::Box) =
    AffineMapping(Diagonal(0.5 * sides(b)), coordinates(center(b)), IHat^pdim(b))

Base.minimum(b::Box) = b.min
Base.maximum(b::Box) = b.max
Base.extrema(b::Box) = b.min, b.max
Base.isapprox(b1::Box, b2::Box) = b1.min ≈ b2.min && b1.max ≈ b2.max
Base.in(p::Point, b::Box) = b.min ≤ p ≤ b.max
