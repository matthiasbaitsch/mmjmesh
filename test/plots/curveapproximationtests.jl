using Test

using MMJMesh
using MMJMesh.Mathematics
using MMJMesh.Plots

import MMJMesh.Plots: PP, Segment, CurveApproximation, refine!, connect!, prev, next, hasprev, hasnext, mark!


# PP tests
p1 = PP(1, 2)
@test p1.x == 1.0
@test p1.point == [1, 2]


# Segment tests
p1 = PP(0, 0)
p2 = PP(1, 0)
p3 = PP(2, 0)
p4 = PP(3, 1)
s1 = Segment(p1, p2)
s2 = Segment(p2, p3)
s3 = Segment(p3, p4)
connect!(s1, s2)
connect!(s2, s3)

@test !hasprev(s1)
@test next(s1) === s2
@test prev(s2) === s1
@test next(s2) === s3
@test prev(s3) === s2
@test !hasnext(s3)
@test angle(prev(s1), s1) == 0
@test angle(s1, next(s1)) == 0
@test angle(s2, next(s2)) ≈ 45
@test angle(s3, next(s3)) == 0

s11 = refine!(s1, PP(0.5, 0.5))
@test s1.level == 1
@test s11.level == 1
@test next(s1) === s11
@test prev(s11) === s1
@test next(s11) === s2
@test prev(s2) === s11

# CurveApproximation tests
f = x -> x^2
p = CurveApproximation(range(-1, 1, 4), f)

@test length(p) == 3
@test mark!(p, 10, 2)
refine!(p, f)
@test length(p) == 6
@test mark!(p, 20, 2)
refine!(p, f)
@test length(p) == 10
@test !mark!(p, 20, 2)

