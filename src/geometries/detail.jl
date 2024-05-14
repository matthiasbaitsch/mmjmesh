mutable struct PointList
    n::Int
    coordinates::Matrix{Float64}
end

Base.iterate(pl::Geometries.PointList) = (pl.coordinates[:, 1], 2)
Base.iterate(pl::Geometries.PointList, i::Int) = i <= pl.n ? (pl.coordinates[:, i], i + 1) : nothing
Base.length(pl::PointList) = pl.n