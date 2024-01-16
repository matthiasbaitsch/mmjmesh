mutable struct ArrayScanner    
    const array::Vector{Any}
    p::Int
    ArrayScanner(array) = new(array, 1)
end

function next!(as::ArrayScanner)
    as.p += 1
    return as.array[as.p - 1]
end

function next!(as::ArrayScanner, n::Int)
    a = as.array[as.p:as.p + n - 1]
    as.p += n
    return a
end

function nextarray!(as::ArrayScanner)
    n = next!(as)
    a = as.array[as.p:as.p + n - 1]
    as.p += n
    return a
end

Base.show(::IO, as::ArrayScanner) = print(io, as.array[as.p:end])
