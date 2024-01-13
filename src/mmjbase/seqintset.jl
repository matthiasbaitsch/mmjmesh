"""
SeqIntSet(a::AbstractVector{Int}; sorted=false)

Set of integers which handles sets containing sequences, for example ```{-100, -99, -98, … , 2, 3, 9, 10, 12, … , 98, 100}```, efficiently.
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

struct Seq
    first::Int
    inc::Int
    last::Int
    Seq(first::Int, inc::Int, last::Int) = new(first, inc, last - (last - first) % inc)
end

Base.first(s::Seq) = s.first
Base.last(s::Seq) = s.last
Base.in(i::Int, s::Seq) = s.first <= i <= s.last && (i - s.first) % s.inc == 0

function Base.show(io::IO, s::Seq)
    if s.first == s.last
        print(io, s.first)
    elseif s.inc == 1
        print(io, s.first, ":", s.last)
    else
        print(io, s.first, ":", s.inc, ":", s.last)
    end
end


nsequences(set::SeqIntSet) = length(set.sequences)
function sequence(set::SeqIntSet, i::Int) 
    s = set.sequences[i]
    return Seq(s.first, 1, s.second)
end

function sequences(s::SeqIntSet)
    return [Seq(p.first, 1, p.second) for p ∈ s.sequences]
end

Base.first(set::SeqIntSet) = first(first(set.sequences))
Base.length(s::SeqIntSet) = s.count[end]
Base.size(s::SeqIntSet) = (length(s),)
Base.isempty(s::SeqIntSet) = (length(s) == 0)
Base.in(target::Int, set::SeqIntSet) = sequenceofvalue(target, set) != -1
Base.getindex(s::SeqIntSet, i::Int) = collect(s)[i]

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, set::SeqIntSet) 
    print(io, "SeqIntSet(ns=$(nsequences(set)), ne=$(length(set)))\n")
    print(io, join(sequences(set), "\n"))
end
Base.show(io::IO, set::SeqIntSet) = print(io, "[$(join(sequences(set), ", "))]")

# TODO: direct access
Base.first(set::SeqIntSet, idx::Int) = first(sequence(set, idx))
inc(set::SeqIntSet, idx::Int) = sequence(set, idx).inc
Base.last(set::SeqIntSet, idx::Int) = last(sequence(set, idx))
Base.last(set::SeqIntSet) = last(set.sequences[end])

# Iterator using the pair (position, value) as state
Base.eltype(::SeqIntSet) = Int
Base.iterate(set::SeqIntSet) = isempty(set) ? nothing : (first(set), (1, first(set)))
function Base.iterate(set::SeqIntSet, state)
    pos, val = state
    if pos <= nsequences(set) && val < last(set)
        if val == last(set, pos)
            pos += inc(set, pos)
            val = first(set, pos)
        else
            val += inc(set, pos)
        end
        return (val, (pos, val))
    else
        return nothing
    end
end

