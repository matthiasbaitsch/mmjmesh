using StaticArrays
using BenchmarkTools

using MMJMesh.Mathematics


# -------------------------------------------------------------------------------------------------
# Monomial evaluation
# -------------------------------------------------------------------------------------------------

function monomialsat1(exponents::SMatrix{N,NT}, x::InRⁿ{N}) where {N,NT}
    c = MVector{NT,Float64}(undef)
    for i = 1:NT
        v = 1
        for j = 1:N
            v *= Base.FastMath.pow_fast(x[j], exponents[j, i])
        end
        c[i] = v
    end
    return c
end

function monomialsat2(exponents::SMatrix{N,NT}, x::InRⁿ{N}) where {N,NT}
    c = Vector{eltype(x)}(undef, NT)
    for i = 1:NT
        v = 1
        for j = 1:N
            v *= Base.FastMath.pow_fast(x[j], exponents[j, i])
        end
        c[i] = v
    end
    return SizedVector{NT}(c)
end

function monomialsat3(exponents::SMatrix{N,NT}, x::InRⁿ{N}) where {N,NT}
    return [
        prod([Base.FastMath.pow_fast(x[j], exponents[j, i]) for j = 1:N])
        for i = 1:NT
    ]
end

e = SA[1 2 3 4 5; 4 3 3 2 1; 1 2 1 2 1]
x = SA[1.2, 2.1, 1.1]

@btime monomialsat1(e, x);
@btime monomialsat2(e, x);
@btime monomialsat3(e, x);



# -------------------------------------------------------------------------------------------------
# MPolynomial evaluation
# -------------------------------------------------------------------------------------------------

# Evaluate monomials
m = _monomialsat(e, x);
println("Evaluation of monomials")
@btime _monomialsat(e, x);

# Scalar valued 
c1 = SA[4, 3, 2, 2, 9]
p1 = MPolynomial(e, c1);
@btime m ⋅ c1;
@btime valueat(p1, x);

# Vector valued 
c2 = SA[4 3 2 2 9; 5 4 5 4 3; 7 7 1 2 9];
p2 = MPolynomial(e, c2);
@btime c2 * m;
@btime valueat(p2, x);