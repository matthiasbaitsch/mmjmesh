using Test
using Random

using MMJMesh
using MMJMesh.Mathematics
using MMJMesh.Plots
import MMJMesh.Plots: sampleadaptive


# -------------------------------------------------------------------------------------------------
# Basic functionality with f(x) = -1 + x^2
# -------------------------------------------------------------------------------------------------

f = Polynomial(-1, 0, 1)
xref = [
    0, 0.07091785989391886, 0.12970579771087964, 0.1729089278540183, 0.23513424487664542, 0.2952937162180036, 0.35851837190925284, 0.43563005931981624, 0.4910605568386503, 0.5680873817855641, 0.636162431409554, 0.6858477591983796, 0.7673958740113812, 0.8132050384999352, 0.8866948711800414, 0.9450280248177776, 1
]

# Number of points
for d ∈ 0:3, np ∈ 3:6
    p = sampleadaptive(f, 0, 1, npoints=np, maxrecursion=d, maxangle=0)
    @test size(p, 2) == (np - 1)^(d + 1) + 1
end

# Compare to reference values
p = sampleadaptive(f, 0, 1, maxrecursion=1, maxangle=0)
@test p[1, :] ≈ xref

# Test scaling factor
p1 = sampleadaptive(f, 0, 1, maxrecursion=20)
p2 = sampleadaptive(41 * f, 0, 1, maxrecursion=20)
p3 = sampleadaptive(41 * f, 0, 1, maxrecursion=20, yscale=1.0 / 41.0)

@test size(p1, 2) < size(p2, 2)
@test p1[1, :] == p3[1, :]

# Insert root
p1 = sampleadaptive(f, -2, 2, maxrecursion=1)
p2 = sampleadaptive(f, -2, 2, maxrecursion=1, ir=true)

@test size(p2, 2) == size(p1, 2) + 2
@test 0 ∈ p2[2, :]
@test p1[:, 1] == p2[:, 1]
@test p1[:, end] == p2[:, end]



