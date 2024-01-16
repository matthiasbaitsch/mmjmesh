# -------------------------------------------------------------------------------------------------
# Geometry struct and constructors
# -------------------------------------------------------------------------------------------------

"""
    Geometry(d)
    Geometry(coordinates)

    Geometry of a finite element mesh.
"""
struct Geometry{D}
    points::PointList
end

Geometry(d::Int, nn::Int=0) = Geometry{d}(PointList(nn, zeros(d, nn)))

Geometry(coordinates::Matrix) =
    Geometry{size(coordinates, 1)}(
        PointList(
            size(coordinates, 2),
            convert(Matrix{Float64}, coordinates)
        ))


# -------------------------------------------------------------------------------------------------
# General methods
# -------------------------------------------------------------------------------------------------

dimension(::Geometry{D}) where {D} = D

"""
    length(g::Geometry, d::int)

Number of geometric objects of dimension `d` in `g`.
"""
function Base.length(g::Geometry{D}, d::Int) where {D}
    @assert 0 <= d <= D
    if d == 0
        return g.points.n
    else
        @notimplemented
    end
end

"""
    g::Geometry[d, idx]

Geometric object of dimension `d` at index `idx` in `g`.
"""
function Base.getindex(g::Geometry{D}, d::Int, idx::Int) where {D}
    @assert 0 <= d <= D
    if d == 0
        return Point{D}(g.points.coordinates[:, idx])
    else
        @notimplemented
    end
end

"""
    squeeze!(g)

Release unused memory from `g`.
"""
function squeeze!(g::Geometry)
    g.coordinates = g.points.coordinates[:, 1:g.points.n]
    return nothing
end


# -------------------------------------------------------------------------------------------------
# Points
# -------------------------------------------------------------------------------------------------

coordinates(g::Geometry) = g.points.coordinates[:, 1:g.points.n]
coordinates(g::Geometry, idx) = g.points.coordinates[:, idx]

function Base.setindex!(g::Geometry{D}, p::AbstractVecOrMat{<:Number}, d::Int, idx::Int) where {D}
    @assert d == 0
    @assert 1 <= idx <= g.points.n
    g.points.coordinates[:, idx] = p
end

Base.setindex!(g::Geometry{D}, p::Point{D}, d::Int, idx::Int) where {D} = setindex!(g, p.coordinates, d, idx)

function Base.push!(g::Geometry{D}, p::Vector{T}) where {D,T<:Number}
    @assert length(p) == D
    if g.points.n == size(g.points.coordinates, 2)
        g.points.coordinates = hcat(g.points.coordinates, Matrix{Float64}(undef, D, g.points.n))
    end
    g.points.n += 1
    g.points.coordinates[:, g.points.n] = p
    return nothing
end

function Base.push!(g::Geometry{D}, pts::Matrix{T}) where {D,T<:Number}
    @assert size(pts, 1) == D
    np = size(pts, 2)
    if g.points.n + np > size(g.points.coordinates, 2)
        nnew = max(g.points.n, np)
        g.points.coordinates = hcat(g.points.coordinates, Matrix{Float64}(undef, D, nnew))
    end
    g.points.coordinates[:, g.points.n+1:g.points.n+np] = pts
    g.points.n += np
    return nothing
end


# -------------------------------------------------------------------------------------------------
# IO
# -------------------------------------------------------------------------------------------------

function Base.show(io::IO, g::Geometry{D}) where D
    println(io, "Geometry{$D} with $(length(g, 0)) Points")
    for i ∈ 1:D
        println(io, "x$i: $(g.points.coordinates[i, 1:g.points.n])")
    end    
end


