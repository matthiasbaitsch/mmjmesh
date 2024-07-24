using Test
using StaticArrays
using LinearAlgebra

using MMJMesh
using MMJMesh.Plots
using MMJMesh.Mathematics

include("validatemappings.jl")


# -------------------------------------------------------------------------------------------------
# On IHat
# -------------------------------------------------------------------------------------------------

C1 = [1, 2]
C2 = [1.0 2.0; 3.0 4.0]
L11 = MappingFromComponents(nodalbasis(makeelement(:lagrange, IHat, k=1))...)

# Interpolate numbers
m1 = Interpolation(L11, C1)
@test m1(-1) == 1.0
@test m1(0) == 1.5
@test m1(1) == 2.0
validate(m1, atol=1e-7)

# Interpolate vectors
m2 = Interpolation(L11, C2)
n = UnitNormal(m2)
@test m2(-1) == [1, 3]
@test m2(0) == [1.5, 3.5]
@test m2(1) == [2, 4]
validate(m2, atol=1e-6)
@test n(0) ≈ [-sqrt(0.5), sqrt(0.5)]


# -------------------------------------------------------------------------------------------------
# On QHat
# -------------------------------------------------------------------------------------------------

C1 = [1, 2, 3, 4]
C2 = [1.0 2.0 2.5 1.1; 0.1 0.2 0.9 1.0]
L21 = MappingFromComponents(nodalbasis(makeelement(:lagrange, QHat, k=1))...)

# Interpolate numbers
m1 = Interpolation(L21, C1)
@test m1(-1, -1) == 1
@test m1(1, -1) == 2
@test m1(1, 1) == 3
@test m1(-1, 1) == 4
@test m1(0, 0) == 2.5
# TODO validate(m1, atol=1e-7)

# Interpolate vectors
m2 = Interpolation(L21, C2)
@test m2(-1, -1) == [1.0, 0.1]
@test m2(1, -1) == [2.0, 0.2]
@test m2(1, 1) == [2.5, 0.9]
@test m2(-1, 1) == [1.1, 1.0]
# TODO validate(m2, atol=1e-6)
