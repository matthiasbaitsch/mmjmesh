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

function _validatederivatives(m::MappingFromR, atol::Real, rtol::Real)
    ed1 = m
    for n ∈ 1:5
        fd = central_fdm(8, n)
        ed1 = ed1'
        ed2 = derivative(m, n)
        for x ∈ _sample(domain(m))
            d1 = derivativeat(m, x, n)
            d2 = ed1(x)
            d3 = ed2(x)
            d4 = fd(m, x)
            @test ed1 == ed2
            @test d1 ≈ d2
            @test d1 ≈ d3
            @test isapprox(d1, d4, atol=atol, rtol=rtol)
        end
    end
end

function _validateantiderivative(f::FunctionRToR)
    fd = central_fdm(8, 1)
    F = antiderivative(f)
    for x ∈ _sample(domain(f))
        v1 = f(x)
        v2 = fd(F, x)
        @test v1 ≈ v2
    end
end

function validate(f::MappingFromR; atol::Real=0.0, rtol::Real=1e-5)

    # Derivatives should always work
    _validatederivatives(f, atol, rtol)

    # Antiderivative should work if it exists
    if hasmethod(antiderivative, Tuple{typeof(f)})
        _validateantiderivative(f)
    end
end

