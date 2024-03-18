function pdim end
function gdim end

atol(Float64) = 1e-10

@enum FromType ROWS COLS
tomatrix(a::Vector{Vector{T}}, from::FromType=COLS) where {T} =
    from == COLS ? stack(a) : stack(a)'

