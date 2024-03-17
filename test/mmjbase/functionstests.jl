using Test

using MMJMesh
using MMJMesh.MMJBase


# -------------------------------------------------------------------------------------------------
# To matrix
# -------------------------------------------------------------------------------------------------

# Input
x = [rand(3) for _ in 1:10]

# From column vectors
a1 = tomatrix(x)
for i ∈ axes(a1, 2)
    @test a1[:, i] == x[i]
end
@test a1 == tomatrix(x, COLS)

# From row vectors
@test tomatrix(x, ROWS) == a1'

