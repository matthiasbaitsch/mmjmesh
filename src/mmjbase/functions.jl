function pdim(o) end
function gdim(o) end

atol(Float64) = 1e-10

@enum FromType ROWS COLS
tomatrix(a::Vector, from::FromType=COLS) =
    from == COLS ? stack(a) : stack(a)'

