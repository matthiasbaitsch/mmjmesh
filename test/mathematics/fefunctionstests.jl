using Test
using StaticArrays
using LinearAlgebra

using MMJMesh
using MMJMesh.Plots
using MMJMesh.Mathematics

include("validatemappings.jl")

# Coefficients
C1 = [1, 2]
C2 = [1.0 2.0; 3.0 4.0]

# Helpers
hx = [3.0, 2.0]
@test MMJMesh.Mathematics.mul(C1, hx) == C1 ⋅ hx
@test MMJMesh.Mathematics.mul(C2, hx) == C2 * hx

# Interpolate numbers
m1 = Interpolation(C1, 1)
@test m1(-1) == 1.0
@test m1(0) == 1.5
@test m1(1) == 2.0
validate(m1, atol=1e-7)

# Interpolate vectors
m2 = Interpolation(C2, 1)
@test m2(-1) == [1, 3]
@test m2(0) == [1.5, 3.5]
@test m2(1) == [2, 4]
validate(m2, atol=1e-6)
