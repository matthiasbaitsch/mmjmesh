"""
SeqIntSet(a::AbstractVector{Int}; sorted=false)

Set of integers which handles sets containing sequences like ```{-100, -99, … , 3, 9, 10, … , 1001}``` efficiently.
"""
struct SeqIntSet
    length::Int
    sequences::Vector{Pair{Int,Int}}
end

function SeqIntSet(a::AbstractVector{Int}; sorted=false)
    isempty(a) && return SeqIntSet(0, Vector{Pair{Int,Int}}())
    asorted = sorted ? a : sort(a)
    entries = Vector{Pair{Int,Int}}()
    start = stop = popfirst!(asorted)
    for i ∈ eachindex(asorted)
        if asorted[i] - stop == 1
            stop = asorted[i]
        else
            push!(entries, Pair(start, stop))
            start = stop = asorted[i]
        end
    end
    push!(entries, Pair(start, stop))
    length = reduce((s, p) -> s + p.second - p.first + 1, entries, init=0)
    return SeqIntSet(length, entries)
end

Base.length(s::SeqIntSet) = s.length
Base.isempty(s::SeqIntSet) = s.length == 0

function Base.in(target::Int, set::SeqIntSet)
    low = 1
    high = length(set.sequences)
    while low <= high
        mid = low + div(high - low, 2)
        seq = set.sequences[mid]
        if seq.first <= target && target <= seq.second
            return true
        elseif seq.second < target
            low = mid + 1
        else
            high = mid - 1
        end
    end
    return false
end

function Base.show(io::IO, s::SeqIntSet)
    if (isempty(s))
        print(io, "[]")
    else
        tostring(p) = p.first == p.second ? string(p.first) : string(p.first) * ":" * string(p.second)
        array = tostring.(s.sequences)
        pstring = "[" * popfirst!(array)
        for a ∈ array
            pstring *= ", " * a
        end
        pstring *= "]"
        print(io, pstring)
    end
end

# Iterator using the pair (position, value) as state
Base.eltype(::SeqIntSet) = Int
Base.iterate(set::SeqIntSet) = isempty(set) ? nothing : (set.sequences[1].first, (1, set.sequences[1].first))
function Base.iterate(set::SeqIntSet, state)
    pos, val = state
    if pos <= length(set.sequences) && val < set.sequences[end].second
        if val == set.sequences[pos].second
            pos += 1
            val = set.sequences[pos].first
        else
            val += 1
        end
        return (val, (pos, val))
    else
        return nothing
    end
end

