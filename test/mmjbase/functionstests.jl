using Test

using MMJMesh
using MMJMesh.MMJBase


# -------------------------------------------------------------------------------------------------
# tomatrix
# -------------------------------------------------------------------------------------------------

vectors = [rand(18) for _ = 1:21]
A = tomatrix(vectors)

for (i, col) = enumerate(eachcol(A))
    @test col == vectors[i]
end
@test tomatrix(vectors, COLS) == A
@test tomatrix(vectors, ROWS) == A'

