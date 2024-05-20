struct ValueAtLF <: AbstractMapping{FunctionRToR,Real,Any}
    x::Real
end
valueat(u::ValueAtLF, f::FunctionToR) = f(u.x)

struct DerivativeAtLF <: AbstractMapping{FunctionRToR,Real,Any}
    x::Real
end
valueat(u::DerivativeAtLF, f::FunctionRToR) = f'(u.x)
