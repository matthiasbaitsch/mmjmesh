
valueat(m::MappingFromRn{N}, x::Vector{<:Real}) where {N} = valueat(m, InRⁿ{N}(x))
(m::MappingFromRn{N})(x::Vector{<:Real}) where {N} = valueat(m, x)
(m::MappingFromRn{N})(x::Real...) where {N} = valueat(m, InRⁿ{N}(x))

derivativeat(m::MappingFromRn{N}, x::Vector{<:Real}, n=1) where {N} =
    derivativeat(m, InRⁿ{N}(x), n)

for func in (:gradientat, :hessianat, :laplacianat)
    @eval begin
        $func(f::FunctionRnToR{N}, x::Vector{<:Real}) where {N} = $func(f, InRⁿ{N}(x))
        $func(f::FunctionRnToR{N}, x::Real...) where {N} = $func(f, InRⁿ{N}(x))
    end
end

divergenceat(v::VectorField{N}, x::Vector{<:Real}) where {N} = divergenceat(v, InRⁿ{N}(x))
divergenceat(v::VectorField{N}, x::Real...) where {N} = divergenceat(v, InRⁿ{N}(x))
