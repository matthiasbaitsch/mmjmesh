using Test

using MMJMesh
using MMJMesh.Mathematics
using MMJMesh.Mathematics: P, Q, S, dimension, _dimension, basis

p3 = Polynomial(1, 2, 3, 4)
p11 = MPolynomial([0 1; 1 0], [1, 2])
p23 = MPolynomial([1 2; 3 2], [1, 2])

@test dimension(Q{1,4}()) == 5
@test dimension(Q{1,4}) == 5
@test dimension(Q{2,4}) == 25
@test dimension(P{2,3}) == 10
@test dimension(P{2,4}) == _dimension(P{2,4})
@test dimension(S{2,2}) == 8
@test dimension(S{2,3}) == 12

@test [1] ∉ Q{2,1}
@test [1, 1] ∈ Q{2,1}
@test [1, 1, 1] ∉ Q{2,1}

@test p3 ∈ Q{1,3}
@test p3 ∈ Q{1,4}
@test p3 ∉ Q{1,2}
@test p3 ∉ Q{2,2}
@test p23 ∉ P{1,10}
@test p23 ∈ P{2,5}
@test p23 ∈ P{2,4}
@test p23 ∉ P{2,3}
@test p23 ∉ P{3,10}
@test p23 ∉ Q{1,10}
@test p23 ∈ Q{2,3}
@test p23 ∈ Q{2,4}
@test p23 ∉ Q{2,2}
@test p23 ∉ Q{3,10}
@test p11 ∈ S{2,2}
@test p23 ∉ S{2,2}
@test p23 ∉ S{2,3}

@test basis(P{1,3}, IHat) == monomials(0:3, IHat)
@test basis(P{1,3}(), IHat) == monomials(0:3, IHat)
@test basis(Q{2,1}, QHat) == mmonomials(2,1,QHat)


