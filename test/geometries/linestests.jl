using Test

using MMJMesh
using MMJMesh.Geometries

# Hline
l = HLine(0.3)
@test parameterof(l, [1.1, 0.3]) == 1.1
@test parameterof(l, [1.1, 0.3 + 1e-13]) == 1.1
@test isnan(parameterof(l, [1.1, 0.4]))
@test [1.1, 0.3] ∈ l
@test [1.1, 0.4] ∉ l

# Vline
l = VLine(0.3)
@test parameterof(l, [0.3, 1.1]) == 1.1
@test parameterof(l, [0.3 + 1e-13, 1.1]) == 1.1
@test isnan(parameterof(l, [0.4, 1.1]))
@test [0.3, 1.1] ∈ l
@test [0.4, 1.1] ∉ l

# Segment
p1 = [3, 2]
p2 = [1, 4]
m = 0.5 * (p1 + p2)
me = m + [1e-12, 1e-12]
s = Segment(p1, p2)

@test parameterof(s, p1) == 0
@test parameterof(s, p2) == 1
@test parameterof(s, m) == 0.5
@test parameterof(s, me, atol=1e-11) == 0.5
@test isnan(parameterof(s, me))
@test isnan(parameterof(s, [5, 5]))
@test m ∈ s
@test [5, 5] ∉ s
@test me ∉ s
@test in(me, s, atol=1e-11)
