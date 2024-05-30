struct ValueAtLF <: AbstractMapping{AbstractMapping,Real,Any}
    x
end
valueat(u::ValueAtLF, f::FunctionToR) = f(u.x)

struct DerivativeAtLF <: AbstractMapping{AbstractMapping,Real,Any}
    x
end
valueat(u::DerivativeAtLF, f::FunctionRToR) = f'(u.x)
