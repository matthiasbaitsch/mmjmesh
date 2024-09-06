using Test
using Symbolics
using IntervalSets
using DomainSets: ×

using MMJMesh
using MMJMesh.Mathematics
using MMJMesh.Mathematics: _combineforms


# -------------------------------------------------------------------------------------------------
# Test _combineforms
# -------------------------------------------------------------------------------------------------

@test _combineforms(points(IHat, :corners), ValueAtLF) |> length == 2
@test _combineforms(
    [points(IHat, :corners), points(IHat, :sides, 2), points(IHat, :interior, 2)],
    ValueAtLF
) |> length == 4
@test _combineforms(points(QHat, :corners), ValueAtLF) |> length == 4
@test _combineforms(points(QHat, :corners), [ValueAtLF, ∂xLF, ∂yLF]) |> length == 12
@test _combineforms(
    [points(QHat, :corners), points(QHat, :sides, 0), points(QHat, :interior, 0)],
    [ValueAtLF, ∂xLF, ∂yLF]
) |> length == 12


# -------------------------------------------------------------------------------------------------
# Test makeelement
# -------------------------------------------------------------------------------------------------

e1 = makeelement(:lagrange, QHat, k=5);
e2 = makeelement(:lagrange, QHat, k=5);
e3 = makeelement(:lagrange, QHat, k=8);
e4 = makeelement(:lagrange, QHat, k=8);
e5 = makeelement(:hermite, QHat);
e6 = makeelement(:hermite, QHat, conforming=false);
@test e1 === e2
@test e1 != e3
@test e3 === e4
@test e5 != e6


# -------------------------------------------------------------------------------------------------
# Validate elements
# -------------------------------------------------------------------------------------------------

function validate(e::FiniteElement; atol)
    N = dimension(e.P)
    ϕ = nodalbasis(e)
    for i = 1:N, j = 1:N
        @test isapprox(e.N[i](ϕ[j]), i == j, atol=atol)
    end
end

validate(makeelement(:lagrange, IHat, k=1), atol=1e-14)
validate(makeelement(:lagrange, QHat, k=1), atol=1e-14)
validate(makeelement(:lagrange, QHat, k=8), atol=1e-8)
validate(makeelement(:serendipity, QHat, k=2), atol=1e-14)
validate(makeelement(:hermite, IHat), atol=1e-14)
validate(makeelement(:hermite, QHat), atol=1e-14)
validate(makeelement(:hermite, QHat, conforming=false), atol=1e-14)


# -------------------------------------------------------------------------------------------------
# Nodal basis
# -------------------------------------------------------------------------------------------------

e = makeelement(:serendipity, QHat, k=2)
ϕ = nodalbasis(e)

# Check that we get the domain right
@test MMJMesh.Mathematics.domain(ϕ[1]) == QHat

# Test nodal property
for i = 1:8, j = 1:8
    @test isequal(e.N[i](ϕ[j]), i == j)
end


# -------------------------------------------------------------------------------------------------
# Handle symbolic domain
# -------------------------------------------------------------------------------------------------

@variables a, b
K = (0 .. a) × (0 .. b)
e = makeelement(:serendipity, K, k=2)
ϕ = nodalbasis(e)

# Check that we get the domain right
@test MMJMesh.Mathematics.domain(ϕ[1]) == R2

# Test nodal property
for i = 1:8, j = 1:8
    @test isequal(e.N[i](ϕ[j]), i == j)
end
