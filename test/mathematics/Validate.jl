module Validate

using Test
using IntervalSets

using MMJMesh
using MMJMesh.Geometries
using MMJMesh.Mathematics

import Random
import DomainSets
import FiniteDifferences

export validate, validateoperations


function _sample(I::Interval)
    xi = [
        0.077555809, 0.095330104, 0.488611758, 0.490830944, 0.541245779,
        0.576334723, 0.648839469, 0.696620560, 0.711463162, 0.782791149
    ]
    a, b = endpoints(I)
    isinf(a) && (a = -5.0)
    isinf(b) && (b = 5.0)
    return a .+ (b - a) * xi
end


function _sample(r::DomainSets.Rectangle)
    @assert isfinite(r)
    Random.seed!(1234)
    N = DomainSets.dimension(r)
    F = parametrization(Box(r))
    return F.([-1 .+ 2 * rand(N) for _ = 1:10])
end


# Validate derivatives for function R to R
function _validatederivatives(m::MappingFromR, atol::Real, rtol::Real)
    ed1 = m
    for n ∈ 1:5
        fd = FiniteDifferences.central_fdm(8, n)
        ed1 = ed1'
        ed2 = derivative(m, n)

        @test ed1 == ed2

        for x ∈ _sample(domain(m))
            d1 = derivativeat(m, x, n)
            d2 = ed1(x)
            d3 = ed2(x)
            d4 = fd(m, x)
            @test d1 ≈ d2
            @test d1 ≈ d3
            @test isapprox(d1, d4, atol=atol, rtol=rtol)
        end
    end
end


# Validate derivatives for functions Rn to R
function _validatederivatives(f::FunctionRnToR, atol::Real, rtol::Real)

    # Validate gradient
    ed1 = f'
    ed2 = gradient(f)
    @test ed1 == ed2
    for x ∈ _sample(domain(f))
        d1 = gradientat(f, x)
        d2 = ed1(x)
        d3 = ed2(x)
        d4 = FiniteDifferences.grad(FiniteDifferences.central_fdm(2, 1), f, Vector(x))[1]
        @test d1 ≈ d2
        @test d1 ≈ d3
        @test isapprox(d1, d4, atol=atol, rtol=rtol)
    end

    # Validate Hessian
    _validatederivatives(gradient(f), atol, rtol)
end


# Validate derivatives for mappings Rn to Rm
function _validatederivatives(m::MappingRnToRm, atol::Real, rtol::Real)
    ed1 = m'
    ed2 = jacobian(m)
    @test ed1 == ed2
    for x ∈ _sample(domain(m))
        d1 = jacobianat(m, x)
        d2 = ed1(x)
        d3 = ed2(x)
        d4 = FiniteDifferences.jacobian(FiniteDifferences.central_fdm(2, 1), m, Vector(x))[1]

        @test derivativeat(m, x) == derivativeat(m, x, 1)

        @test d1 ≈ d2
        @test d1 ≈ d3
        @test isapprox(d1, d4, atol=atol, rtol=rtol)
    end
end


function _validateantiderivative(f::MappingFromR)
    fd = FiniteDifferences.central_fdm(8, 1)
    F = antiderivative(f)
    for x ∈ _sample(domain(f))
        v1 = f(x)
        v2 = fd(F, x)
        @test v1 ≈ v2
    end
end

_validatecodomaintype(x, ::Type{<:InR}) = x isa Real
_validatecodomaintype(x, ::Type{<:InRⁿ{N}}) where {N} = (size(x) == (N,))
_validatecodomaintype(x, ::Type{<:InRⁿˣᵐ{N,M}}) where {N,M} = size(x) == (N, M)


function _validatecodomaintype(m::AbstractMapping)
    for x ∈ _sample(domain(m))
        @test _validatecodomaintype(valueat(m, x), codomaintype(m))
    end
    return true
end


function validate(m::AbstractMapping; atol::Real=0.0, rtol::Real=1e-5)

    # Types
    _validatecodomaintype(m)

    # Derivatives
    _validatederivatives(m, atol, rtol)

    # Antiderivative
    if m isa FunctionRToR
        try
            _validateantiderivative(m)
        catch e
            if !(e isa MMJMesh.MMJBase.NotImplementedError)
                rethrow(e)
            end
        end
    end
    return true
end


function validate(e::FiniteElement; atol)
    N = dimension(e.P)
    ϕ = nodalbasis(e)
    for i = 1:N, j = 1:N
        @test isapprox(e.N[i](ϕ[j]), i == j, atol=atol)
    end
end


function validateoperations(f1, f2)
    h1 = f1 + f2
    h2 = f1 * f2
    h3 = π * f1
    for x ∈ _sample(domain(f1) ∩ domain(f2))
        @test h1(x) ≈ f1(x) + f2(x)
        @test h2(x) ≈ f1(x) * f2(x)
        @test h3(x) ≈ π * f1(x)
    end
    return true
end

end
