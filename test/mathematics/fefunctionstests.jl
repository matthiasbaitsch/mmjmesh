using Test
using StaticArrays
using LinearAlgebra

using MMJMesh
using MMJMesh.Plots
using MMJMesh.Mathematics

# Coefficients
C1 = [1, 2]
C2 = [1.0 2.0; 3.0 4.0]

# Interpolate numbers
m1 = Interpolation(C1, 1)
@test m1(-1) == 1.0
@test m1(0) == 1.5
@test m1(1) == 2.0
validate(m1, atol=1e-7)

# Interpolate vectors
m2 = Interpolation(C2, 1)
n = UnitNormal(m2)
@test m2(-1) == [1, 3]
@test m2(0) == [1.5, 3.5]
@test m2(1) == [2, 4]
validate(m2, atol=1e-6)
@test n(0) ≈ [-sqrt(0.5), sqrt(0.5)]
