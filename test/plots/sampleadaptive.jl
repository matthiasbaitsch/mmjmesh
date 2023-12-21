using Test
using LinearAlgebra

import MMJMesh.Plots: valuerange, approximationerror, safeeval, X1, X2, W1, W2

# TODO: Tests do not make that much sense

# Numerical integration formulae
function nint(f, a, b)
    w = b - a
    y1 = f.(a .+ w * X1)
    y2 = f.(a .+ w * X2)
    w * dot(W1, y1), w * dot(W2, y2)
end
check(a, b) = a[1] ≈ b[1] && a[2] ≈ b[2]
@test check(nint(x -> x, 0, 2), (2.0, 2.0))

# Helper functions
@test isnan(safeeval(sin, 1.0 / 0.0))
@test valuerange(sin, 0, 1, 21) - 2 < 1e-4
@test abs(valuerange(x -> sin(1/x), 0, 1, 21) - 2) < 0.5

# Approximation error
function ae(f, a, b)
    X = X1
    W = W1
    h = 1
    w = b - a
    x = a .+ w * X
    y = f.(x)
    return approximationerror(h, X, W, y)
end
@test abs(ae(sin, 0, 1) - 0.043) < 1e-3
