"""
SeqIntSet(a::AbstractVector{Int}; sorted=false)

Set of integers which handles sets containing sequences, for example ```{-100, -99, -98, … , 2, 3, 9, 10, 12, … , 98, 100}```, efficiently.
"""
struct SeqIntSet{T<:Integer} <: AbstractVector{T}
    start::Vector{T}
    step::Vector{T}
    index::Vector{T}
end

# Constructor
function SeqIntSet(a::AbstractVector{T}) where {T<:Integer}

    # Quick return
    isempty(a) && return SeqIntSet(Int[], Int[], [1])

    # Initialize
    start = T[]
    step = T[]
    index = T[]
    state = 0
    step1 = step2 = NaN
    asu = a |> sort |> unique

    # Read sequence
    for i = 2:length(asu)

        # Calculate steps
        step1 = step2
        step2 = asu[i] - asu[i-1]

        # Process states
        if state == 0
            push!(index, i - 1)
            push!(start, asu[i-1])
            state = 1
        elseif state == 1
            if step1 != step2
                push!(step, step1)
                push!(index, i - 1)
                push!(start, asu[i-1])
            else
                state = 2
            end
        else
            if step1 != step2
                push!(step, step1)
                state = 0
            end
        end
    end

    # Wrap up
    if state <= 1
        if state == 1
            push!(step, step2)
        end
        push!(index, length(asu))
        push!(start, asu[end])
        push!(step, 0)
    else
        push!(step, step1)
    end
    push!(index, length(asu) + 1)

    # Return
    return SeqIntSet(start, step, index)
end

# Helpers
nsequences(set::SeqIntSet) = length(set.start)
nsteps(set::SeqIntSet, s::Int) = set.index[s+1] - set.index[s] - 1
start(set::SeqIntSet, s::Int) = set.start[s]
step(set::SeqIntSet, s::Int) = set.step[s]
stop(set::SeqIntSet, s::Int) = set.start[s] + nsteps(set, s) * set.step[s]
index(set::SeqIntSet, s::Int) = set.index[s]

# Base
Base.first(set::SeqIntSet) = first(set.start)
Base.length(set::SeqIntSet) = set.index[end] - 1
Base.size(set::SeqIntSet) = (length(set),)
Base.isempty(set::SeqIntSet) = (length(set) == 0)
Base.last(set::SeqIntSet) = set.start[end] + nsteps(set, nsequences(set)) * set.step[end]
Base.union(s1::SeqIntSet, s2::SeqIntSet) = SeqIntSet([collect(s1); collect(s2)])
Base.intersect(s1::SeqIntSet, s2::SeqIntSet) = SeqIntSet([i for i in s1 if i ∈ s2])
Base.setdiff(s1::SeqIntSet, s2::SeqIntSet) = SeqIntSet([i for i in s1 if i ∉ s2])

# Getindex
function Base.getindex(set::SeqIntSet, i::Int)
    @assert 1 <= i <= length(set)
    idx = searchsortedlast(set.index, i)
    return set.start[idx] + (i - set.index[idx]) * set.step[idx]
end

# In
function Base.in(target::Int, set::SeqIntSet)
    s = searchsortedlast(set.start, target)
    s == 0 && return false
    target > stop(set, s) && return false
    step(set, s) == 0 && return true
    return (target - start(set, s)) % step(set, s) == 0
end

# Iterate
Base.eltype(::SeqIntSet) = Int
Base.iterate(set::SeqIntSet) = isempty(set) ? nothing : (first(set), (1, 2))
function Base.iterate(set::SeqIntSet, state)
    sequence, index = state
    if sequence <= nsequences(set) && index <= length(set)
        if index == set.index[sequence+1]
            sequence += 1
        end
        offset = index - set.index[sequence]
        value = set.start[sequence] + offset * set.step[sequence]
        return (value, (sequence, index + 1))
    end
    return nothing
end

# Show
_ts(start, step, stop) = start == stop ? string(start) : step == 1 ? string(start, ":", stop) : string(start, ":", step, ":", stop)
_tsseq(set::SeqIntSet) = [_ts(start(set, i), step(set, i), stop(set, i)) for i in 1:nsequences(set)]
function Base.show(io::IO, ::MIME{Symbol("text/plain")}, set::SeqIntSet)
    print(io, "SeqIntSet(ns=$(nsequences(set)), ne=$(length(set)))\n")
    print(io, join(_tsseq(set), "\n"))
end
Base.show(io::IO, set::SeqIntSet) = print(io, "[$(join(_tsseq(set), ", "))]")