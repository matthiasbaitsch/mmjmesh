"""
    ConnectivityList

A `ConnectivityList` stores connections between mesh entities. Each entry in a connectivity is an array
of links to other entities.

TODO: Complete implementation of `AbstractVector` interface.
"""
struct ConnectivityList <: AbstractVector{Vector{Int}}
    offsets::Vector{Int}
    entries::Vector{Int}
end

ConnectivityList() = ConnectivityList([1], Int[])

function ConnectivityList(linkslist::Vector{Vector{Int}})
    cl = ConnectivityList()
    for links in linkslist
        push!(cl, links)
    end
    return cl
end


# -------------------------------------------------------------------------------------------------
# Methods
# -------------------------------------------------------------------------------------------------

function inverse(cl::ConnectivityList)
    m = Dict{Set{Int},Int}()
    for (i, c) in enumerate(cl)
        m[Set(c)] = i
    end
    return m
end

maxlinksize(cl::ConnectivityList) = maximum([length(links) for links in cl])


# -------------------------------------------------------------------------------------------------
# Adapt methods from Base
# -------------------------------------------------------------------------------------------------

Base.:(==)(cl1::ConnectivityList, cl2::ConnectivityList) = cl1.offsets == cl2.offsets && cl1.entries == cl2.entries
Base.length(cl::ConnectivityList) = length(cl.offsets) - 1
Base.length(cl::ConnectivityList, i::Int) = cl.offsets[i+1] - cl.offsets[i]
Base.iterate(cl::ConnectivityList, state=1) = state > length(cl) ? nothing : (cl[state], state + 1)
Base.size(cl::ConnectivityList) = (length(cl), )

function Base.push!(connectivityList::ConnectivityList, links::Vector{Int})
    push!(connectivityList.offsets, connectivityList.offsets[end] + length(links))
    push!(connectivityList.entries, links...)
    return nothing
end

Base.getindex(cl::ConnectivityList, i::Int) = cl.entries[(cl.offsets[i]):(cl.offsets[i+1]-1)]

function Base.getindex(cl::ConnectivityList, i::Int, j::Int)
    @assert 1 <= i <= length(cl) && 1 <= j <= length(cl, i)
    return cl.entries[cl.offsets[i]+j-1]
end

function Base.adjoint(cl::ConnectivityList)
    n = maximum(cl.entries)
    linkslist = [Vector{Int}() for _ in 1:n]
    for (i, links) in enumerate(cl)
        for j in links
            push!(linkslist[j], i)
        end
    end
    return ConnectivityList(linkslist)
end


# -------------------------------------------------------------------------------------------------
# IO
# -------------------------------------------------------------------------------------------------

# function Base.show(io::IO, cl::ConnectivityList)
#     n = length(cl)
#     mnl = maxlinksize(cl)
#     for (i, links) in enumerate(cl)
#         s1 = @sprintf "%7i:" i
#         s2 = join([@sprintf "%8i" l for l in links])
#         if n <= 10 || i <= 5 || i >= n - 5
#             println(io, s1, s2)
#         end
#         if n > 10 && i == 5
#             println(io, "      ⋮ ", repeat("       ⋮", mnl))
#         end
#     end
# end
