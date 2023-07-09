import Base.adjoint
import Base.getindex
import Base.iterate
import Base.length
import Base.push!
import Base.show
import Base.==

using Printf

struct ConnectivityList
    offsets::Vector{Int}
    entries::Vector{Int}
end

function ConnectivityList()
    return ConnectivityList([1], Int[])
end

function ConnectivityList(linksList::Vector{Vector{Int}})
    cl = ConnectivityList()
    for links in linksList
        push!(cl, links)
    end
    return cl
end

(==)(cl1::ConnectivityList, cl2::ConnectivityList) = cl1.offsets == cl2.offsets && cl1.entries == cl2.entries
Base.length(cl::ConnectivityList) = length(cl.offsets) - 1
Base.length(cl::ConnectivityList, i::Int) = cl.offsets[i+1] - cl.offsets[i]
Base.getindex(cl::ConnectivityList, i::Int) = cl.entries[(cl.offsets[i]):(cl.offsets[i + 1] - 1)]
Base.iterate(cl::ConnectivityList, state=1) = state > length(cl) ? nothing : (cl[state], state + 1)

function push!(connectivityList::ConnectivityList, links::Vector{Int})
    push!(connectivityList.offsets, connectivityList.offsets[end] + length(links))
    push!(connectivityList.entries, links...)
    return nothing
end

function Base.getindex(cl::ConnectivityList, i::Int, j::Int)
    @assert 1 <= i <= length(cl)
    @assert 1 <= j <= length(cl, i)
    return cl.entries[cl.offsets[i]+j-1]
end

function Base.adjoint(cl::ConnectivityList)
    n = maximum(cl.entries)
    linksList = [Vector{Int}() for _ in 1:n]
    for (i, links) in enumerate(cl)
        for j in links
            push!(linksList[j], i)
        end
    end
    return ConnectivityList(linksList)
end

function inverse(cl::ConnectivityList) 
    m = Dict{Set{Int}, Int}();
    for (i, c) in enumerate(cl)
        m[Set(c)] = i
    end
    return m
end

function Base.show(io::IO, cl::ConnectivityList)
    for (i, links) in enumerate(cl)
        s1 = @sprintf "%7i:" i
        s2 = join([@sprintf "%10i" l for l in links])
        println(io, s1, s2)
    end
end
