atol(Float64) = 1e-10


# Shorthand for collecting vectors into a matrix
@enum FromType ROWS COLS
tomatrix(a, from::FromType=COLS) = from == COLS ? stack(a) : stack(a, dims=1)


# Unit vector
unitvector(idx::Integer, n::Integer) = [idx == i for i = 1:n]

