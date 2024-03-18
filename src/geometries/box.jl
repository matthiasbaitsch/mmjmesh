struct Box{D} <: GeometricObject{D,D}
    min::Point{D}
    max::Point{D}

    function Box{D}(min, max) where {D}
        @assert min ≤ max "`min` must be less than or equal to `max`"
        new(min, max)
    end
end

Box(min::Point{D}, max::Point{D}) where {D} = Box{D}(min, max)
Box(min::Tuple, max::Tuple) = Box(Point(min...), Point(max...))

gdim(::Type{<:Box{D}}) where {D} = D
pdim(::Type{<:Box{D}}) where {D} = D

center(b::Box) = Point((coordinates(b.max) + coordinates(b.min)) / 2)
diagonal(b::Box) = norm(b.max - b.min)
sides(b::Box) = Tuple(b.max - b.min)

Base.minimum(b::Box) = b.min
Base.maximum(b::Box) = b.max
Base.extrema(b::Box) = b.min, b.max
Base.isapprox(b1::Box, b2::Box) = b1.min ≈ b2.min && b1.max ≈ b2.max
Base.in(p::Point, b::Box) = b.min ≤ p ≤ b.max
