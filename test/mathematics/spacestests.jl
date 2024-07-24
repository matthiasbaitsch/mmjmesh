using Test
using Symbolics
using IntervalSets
using DomainSets: ×

using MMJMesh
using MMJMesh.Mathematics
using MMJMesh.Mathematics: P, Q, S, dimension, _dimension, basis


# -------------------------------------------------------------------------------------------------
# Example spaces
# -------------------------------------------------------------------------------------------------

const P22 = P{2,2,QHat}
const P23 = P{2,3,QHat}
const P24 = P{2,4,QHat}
const P25 = P{2,5,QHat}
const Q11 = Q{1,1,IHat}
const Q12 = Q{1,2,IHat}
const Q13 = Q{1,3,IHat}
const Q14 = Q{1,4,IHat}
const Q21 = Q{2,1,QHat}
const Q22 = Q{2,2,QHat}
const Q23 = Q{2,3,QHat}
const Q24 = Q{2,4,QHat}
const S22 = S{2,2,QHat}
const S23 = S{2,3,QHat}


# -------------------------------------------------------------------------------------------------
# Dimension
# -------------------------------------------------------------------------------------------------

@test dimension(Q{1,4,IHat}()) == 5
@test dimension(Q14) == 5
@test dimension(Q24) == 25
@test dimension(P23) == 10
@test dimension(P24) == _dimension(P24)
@test dimension(S22) == 8
@test dimension(S23) == 12
@test dimension(Q23R) == 12
@test _dimension(Q23R) == 12


# -------------------------------------------------------------------------------------------------
# Function in
# -------------------------------------------------------------------------------------------------

p3 = Polynomial(1, 2, 3, 4)
p11 = MPolynomial([0 1; 1 0], [1, 2])
p23 = MPolynomial([1 2; 3 2], [1, 2])

@test [1] ∉ Q21
@test [1, 1] ∈ Q21
@test [1, 1, 1] ∉ Q21
@test p3 ∈ Q13
@test p3 ∈ Q14
@test p3 ∉ Q12
@test p3 ∉ Q22
@test p23 ∉ P{1,10,IHat}
@test p23 ∈ P{2,5,QHat}
@test p23 ∈ P24
@test p23 ∉ P23
@test p23 ∉ P{3,10,QHat}
@test p23 ∉ Q{1,10,QHat}
@test p23 ∈ Q{2,3,QHat}
@test p23 ∈ Q24
@test p23 ∉ Q22
@test p23 ∉ Q{3,10,QHat}
@test p11 ∈ S22
@test p23 ∉ S22
@test p23 ∉ S23



# -------------------------------------------------------------------------------------------------
# Basis
# -------------------------------------------------------------------------------------------------

@test basis(P{1,3,IHat}) == monomials(0:3, IHat)
@test basis(P{1,3,IHat}()) == monomials(0:3, IHat)
@test basis(Q21) == mmonomials(2, 1, QHat)
@test length(basis(Q23R{QHat})) == 12


# -------------------------------------------------------------------------------------------------
# Symbolic domains
# -------------------------------------------------------------------------------------------------

# TODO - Remove domain???
# @variables a, b
# K = (0 .. a) × (0 .. b)
# QX = Q{2,1,K}()

