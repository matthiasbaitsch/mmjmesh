# -------------------------------------------------------------------------------------------------
# Struct and constructors
# -------------------------------------------------------------------------------------------------

"""
    ConnectivityList

A `ConnectivityList` stores connections between mesh entities. Each entry in a `ConnectivityList` is an array
of links to other entities.
"""
struct ConnectivityList <: AbstractVector{Vector{Int}}
    offsets::Vector{Int}
    entries::Vector{Int}
end

"""
    ConnectivityList()

Construct empty `ConnectivityList`.
"""
ConnectivityList() = ConnectivityList([1], Int[])

"""
    ConnectivityList(linkslist::Vector{Vector{Int}})

Construct `ConnectivityList` from a list of links.
"""
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

"""
    inverse(cl::ConnectivityList)

The inverse of a connectivity list maps links to an index. Note that the inverse does not
take into account the order of the links.
"""
function inverse(cl::ConnectivityList)
    m = Dict{Set{Int},Int}()
    for (i, c) in enumerate(cl)
        m[Set(c)] = i
    end
    return m
end

"""
    maxlinkssize(cl::ConnectivityList)

Get the maximum link size of a connectivity list.
"""
maxlinkssize(cl::ConnectivityList) = maximum([length(links) for links in cl])


# -------------------------------------------------------------------------------------------------
# Methods from Base
# -------------------------------------------------------------------------------------------------

"""
    push!(cl::ConnectivityList, links::Vector{Int})

Add entity by adding links to the end of the list.
"""
function Base.push!(cl::ConnectivityList, links::Vector{Int})
    push!(cl.offsets, cl.offsets[end] + length(links))
    push!(cl.entries, links...)
    return nothing
end

function Base.append!(cl::ConnectivityList, links::AbstractVector{Vector{Int}})
    for l in links
        push!(cl, l)
    end
    return nothing
end

"""
    length(cl)

Returns the number of entities.
"""
Base.length(cl::ConnectivityList) = length(cl.offsets) - 1

"""
    length(cl, i)

Returns the number of links from entity `i`.
"""
Base.length(cl::ConnectivityList, i::Int) = cl.offsets[i+1] - cl.offsets[i]

"""
    cl[i]

Returns the links from entity `i`.
"""
Base.getindex(cl::ConnectivityList, i::Int) = cl.entries[(cl.offsets[i]):(cl.offsets[i+1]-1)]

"""
    cl[i, j]

Returns the `j`-th link from entity `i`.
"""
function Base.getindex(cl::ConnectivityList, i::Int, j::Int)
    @assert 1 <= i <= length(cl) && 1 <= j <= length(cl, i)
    return cl.entries[cl.offsets[i]+j-1]
end

"""
    cl'

Returns the transpose of connectivity list `cl`.
"""
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

Base.:(==)(cl1::ConnectivityList, cl2::ConnectivityList) = cl1.offsets == cl2.offsets && cl1.entries == cl2.entries
Base.size(cl::ConnectivityList) = (length(cl),)
Base.iterate(cl::ConnectivityList, state=1) = state > length(cl) ? nothing : (cl[state], state + 1)


# -------------------------------------------------------------------------------------------------
# IO
# -------------------------------------------------------------------------------------------------

# function Base.show(io::IO, cl::ConnectivityList)
#     n = length(cl)
#     mnl = maxlinkssize(cl)
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

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, cl::ConnectivityList)
    println(io, "$(length(cl))-element ConnectivityList:")
    print(io, " " * join(collect(cl), "\n "))
end

function Base.show(io::IO, cl::ConnectivityList)
    print(io, "[$(join(collect(cl), ", "))]")
end
