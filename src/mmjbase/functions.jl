
atol(Float64) = 1e-10

@enum FromType ROWS COLS
tomatrix(a, from::FromType=COLS) = from == COLS ? stack(a) : stack(a, dims=1)

