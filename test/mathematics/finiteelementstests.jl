@testitem "hat functions" begin

    using Test
    using LinearAlgebra: ⋅

    using MMJMesh
    using MMJMesh.Mathematics

    n = 5
    xs = range(0 .. 4.23, n)
    uhats = 1:n
    phis = hatfunctions(xs)

    uh1 = uhats ⋅ phis
    uh2 = phis ⋅ uhats
    @test uh1.(xs) ≈ uhats
    @test uh2.(xs) ≈ uhats
    @test iszero(zeros(n) ⋅ phis)

end


@testitem "combine forms" begin

    using Test

    using MMJMesh
    using MMJMesh.Mathematics
    using MMJMesh.Mathematics: _combineforms

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

end


@testitem "make element" begin

    using Test

    using MMJMesh
    using MMJMesh.Mathematics

    e1 = makeelement(:lagrange, QHat, k=5)
    e2 = makeelement(:lagrange, QHat, k=5)
    e3 = makeelement(:lagrange, QHat, k=8)
    e4 = makeelement(:lagrange, QHat, k=8)
    e5 = makeelement(:hermite, QHat)
    e6 = makeelement(:hermite, QHat, conforming=false)
    @test e1 === e2
    @test e1 != e3
    @test e3 === e4
    @test e5 != e6

end


@testitem "validate elements" setup = [Validate] begin

    using Test

    using MMJMesh.Mathematics

    Validate.validate(makeelement(:lagrange, IHat, k=1), atol=1e-14)
    Validate.validate(makeelement(:lagrange, QHat, k=1), atol=1e-14)
    Validate.validate(makeelement(:lagrange, QHat, k=8), atol=1e-8)
    Validate.validate(makeelement(:serendipity, QHat, k=2), atol=1e-14)
    Validate.validate(makeelement(:hermite, IHat), atol=1e-14)
    Validate.validate(makeelement(:hermite, QHat), atol=1e-14)
    Validate.validate(makeelement(:hermite, QHat, conforming=false), atol=1e-14)

end


@testitem "nodal basis" begin

    using Test

    using MMJMesh
    using MMJMesh.Mathematics

    e = makeelement(:serendipity, QHat, k=2)
    ϕ = nodalbasis(e)

    # Check that we get the domain right
    @test MMJMesh.Mathematics.domain(ϕ[1]) == QHat

    # Test nodal property
    for i = 1:8, j = 1:8
        @test isequal(e.N[i](ϕ[j]), i == j)
    end

    @test nodalbasis(:serendipity, QHat; k=2) == ϕ

end


@testitem "symbolic domain" begin

    using Test
    using Symbolics

    using MMJMesh
    using MMJMesh.Mathematics
    using DomainSets: ×

    @variables a, b
    K = (0 .. a) × (0 .. b)

    # Serendipity
    e = makeelement(:serendipity, K, k=2)
    ϕ = nodalbasis(e)

    @test MMJMesh.Mathematics.domain(ϕ[1]) == R2
    for i = 1:8, j = 1:8
        @test isequal(e.N[i](ϕ[j]), i == j)
    end

end
