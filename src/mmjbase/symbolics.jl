function Base.rationalize(expression::Symbolics.Num; tol::Real=eps(Float32))
    Symbolics.@variables zero, one
    dorationalize(x) = false
    dorationalize(x::AbstractFloat) = true
    rule = Symbolics.@rule ~x::dorationalize => (rationalize(~x, tol=tol) + zero)
    rewriter = SymbolicUtils.Postwalk(Symbolics.Chain([rule]))
    expression = Symbolics.simplify(expression, rewriter=rewriter)
    expression = Symbolics.simplify(Symbolics.substitute(expression, Dict(zero => 0, one => 1)))
    return expression
end

rationalize!(c::AbstractArray{Symbolics.Num}) = map!(rationalize, c, c)

function integerize(expression::Symbolics.Num)
    Symbolics.@variables xone, xnull
    dointegerize(x::Rational) = (denominator(x) == 1)
    dointegerize(x::AbstractFloat) = (x == round(x, digits=0))
    dointegerize(x) = false
    r = Symbolics.@rule ~x::dointegerize => (Int(~x) + xnull)
    expression = Symbolics.simplify(expression)
    expression = Symbolics.simplify(expression, rewriter=SymbolicUtils.Postwalk(Symbolics.Chain([r])))
    expression = Symbolics.simplify(Symbolics.substitute(expression, Dict(xone => 1, xnull => 0)))
    return expression
end

integerize!(c::AbstractArray{Symbolics.Num}) = map!(integerize, c, c)

"""
    integerize(x) -> x

Fallback for non-symbolic values (e.g. a plain `Float64` coefficient): returned unchanged.
"""
integerize(x) = x

"""
    cancel(expr::Symbolics.Num) -> Symbolics.Num

Divide the numerator and denominator of the rational expression `expr` by the greatest
common divisor of all their integer/rational coefficients.

Repeated symbolic fraction arithmetic (as done internally by `integrate`) tends to leave
numerator and denominator with a huge but cancellable common factor that Symbolics' own
`simplify`/`expand` do not remove, since they cancel shared symbolic factors but do not
reduce the numeric content of a multi-term rational expression.
"""
function cancel(expr::Symbolics.Num)

    function _numericcoefficients(expr::Symbolics.Num)
        vars = Symbolics.get_variables(expr)
        isempty(vars) && return Any[Symbolics.value(expr)]
        cs, rem = Symbolics.polynomial_coeffs(Symbolics.value(expr), vars)
        vals = Any[SymbolicUtils.unwrap_const(v) for v in values(cs)]
        remv = SymbolicUtils.unwrap_const(rem)
        remv isa Union{Integer,Rational} && !iszero(remv) && push!(vals, remv)
        return vals
    end


    n, d = numerator(expr), denominator(expr)
    isequal(d, 1) && return expr

    nums = filter(
        v -> v isa Union{Integer,Rational},
        vcat(_numericcoefficients(n), _numericcoefficients(d))
    )
    isempty(nums) && return expr
    g = reduce(gcd, nums)
    (g == 0 || g == 1) && return expr

    return integerize(Symbolics.simplify(n / g; expand=true)) / integerize(Symbolics.simplify(d / g; expand=true))
end

"""
    cancel(x) -> x

Fallback for non-symbolic results (e.g. a plain `Float64` from numeric integration):
returned unchanged.
"""
cancel(x) = x