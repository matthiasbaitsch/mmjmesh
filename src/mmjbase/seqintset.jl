"""
SeqIntSet(a::AbstractVector{Int}; sorted=false)

Set of integers which handles sets containing sequences like ```{-100, -99, … , 3, 9, 10, … , 1001}``` efficiently.
"""
struct SeqIntSet <: AbstractVector{Int}
    count::Vector{Int}
    sequences::Vector{Pair{Int,Int}}
end

function SeqIntSet(a::AbstractVector{Int}; sorted=false)
    isempty(a) && return SeqIntSet([0], Vector{Pair{Int,Int}}())
    asorted = sorted ? a : sort(a)
    sequences = Vector{Pair{Int,Int}}()
    start = stop = popfirst!(asorted)
    for i ∈ eachindex(asorted)
        if asorted[i] - stop <= 1
            stop = asorted[i]
        else
            push!(sequences, Pair(start, stop))
            start = stop = asorted[i]
        end
    end
    push!(sequences, Pair(start, stop))
    count = [p.second - p.first + 1 for p ∈ sequences] |> cumsum
    return SeqIntSet(count, sequences)
end

function sequenceofvalue(target::Int, set::SeqIntSet)
    low = 1
    high = length(set.sequences)
    while low <= high
        mid = low + div(high - low, 2)
        seq = set.sequences[mid]
        if seq.first <= target && target <= seq.second
            return mid
        elseif seq.second < target
            low = mid + 1
        else
            high = mid - 1
        end
    end
    return -1
end

Base.length(s::SeqIntSet) = s.count[end]
Base.size(s::SeqIntSet) = (length(s),)
Base.isempty(s::SeqIntSet) = length(s) == 0
Base.in(target::Int, set::SeqIntSet) = sequenceofvalue(target, set) != -1

# TODO: Implement this reasonably
Base.getindex(s::SeqIntSet, i::Int) = collect(s)[i]

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

