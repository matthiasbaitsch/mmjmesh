abstract type Region{D} end

function in end
function discretize end

"""
    dim(R::Region)

Return the dimension of the region `R`.
"""
dim(::Region{D}) where {D} = D


struct Interval <: Region{1}
    a::Float64
    b::Float64

    function Interval(a, b)
        @assert a <= b
        return new(a, b)
    end
end

Base.isequal(I1::Interval, I2::Interval) = I1.a == I2.a && I1.b == I2.b
Base.isapprox(I1::Interval, I2::Interval; atol=atol(Float64), kwargs...) = isapprox(I1.a, I2.a; atol, kwargs...) && isapprox(I1.b, I2.b; atol, kwargs...)

Base.in(I::Interval, x::Number; eps::Float64 = 0.0) = I.a - eps <= x <= I.b + eps
discretize(I::Interval, n::Int) = [I.a + (I.b - I.a) * i / (n - 1) for i in 0:n-1]

