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