using Test

using MMJMesh
import MMJMesh.MMJBase: SeqIntSet, sequenceofvalue

# Many elements
s = SeqIntSet([1, 10, -50, 2, 11, 9, 8, 0, 100]);

@test sequenceofvalue(-100, s) == -1
@test sequenceofvalue(-51, s) == -1
@test sequenceofvalue(-50, s) == 1
@test sequenceofvalue(0, s) == 2
@test sequenceofvalue(1, s) == 2
@test sequenceofvalue(2, s) == 2
@test sequenceofvalue(8, s) == 3
@test sequenceofvalue(11, s) == 3
@test sequenceofvalue(100, s) == 4
@test sequenceofvalue(101, s) == -1
@test sequenceofvalue(1000, s) == -1

@test length(s.sequences) == 4
@test size(s.sequences) == (4,)
@test length(s) == 9
@test collect(s) == [-50, 0, 1, 2, 8, 9, 10, 11, 100]
@test -1000 ∉ s
@test -51 ∉ s
@test 5 ∉ s
@test 12 ∉ s
@test 101 ∉ s
@test 1000 ∉ s
@test all([i ∈ s for i ∈ s])

# TODO: Efficient getindex
#function Base.getindex(s::SeqIntSet, i::Int)
#end
# function sequenceofindex(idx::Int, set::SeqIntSet)
#     i1 = searchsorted(set.count, idx)

#     i1
# end
# sequenceofindex(1, s)
# sequenceofindex(2, s)
# sequenceofindex(3, s)
# @test sequenceofindex(1, s) == 1
# @test sequenceofindex(2, s) == 2
# @test sequenceofindex(3, s) == 2
# @test sequenceofindex(4, s) == 2
# @test sequenceofindex(5, s) == 3
# @test sequenceofindex(8, s) == 3
# @test sequenceofindex(9, s) == 4

# One element
s = SeqIntSet([1])
@test length(s) == 1
@test collect(s) == [1]
@test 1 ∈ s
@test -1 ∉ s
@test 2 ∉ s

# Empty
s = SeqIntSet(Vector{Int}())
@test isempty(s)
@test collect(s) == Vector{Int}()
@test 5 ∉ s

# Array indexing
a = collect(1:10)
s = SeqIntSet([1, 2, 3, 5, 9, 10])
@test a[s] == s

# TODO: Set operations
# s1 = SeqIntSet([1, 2, 3, 5, 9, 10])
# s2 = SeqIntSet([1, 2, 3, 4, 5, 6])

# Base.union(s1::SeqIntSet, s2::SeqIntSet) = SeqIntSet(union(s1.sequences, s2.sequences))

# @test s1 ∪ s2 == [1, 2, 3, 4, 5, 6, 9, 10]
# @test s1 ∩ s2 == [1, 2, 3, 5]



# Test case from a mesh
a1 = [1, 2, 1, 12, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 22, 12, 23, 22, 33, 23, 34, 33, 44, 34, 45, 44, 55, 45, 56, 55, 66, 56, 67, 66, 77, 67, 78, 77, 88, 78, 89, 88, 99, 89, 100, 99, 110, 100, 111, 110, 121, 111, 122, 121, 132, 122, 133, 132, 143, 133, 144, 143, 154, 144, 155, 154, 165, 166, 167, 155, 166, 167, 168, 168, 169, 169, 170, 170, 171, 171, 172, 172, 173, 173, 174, 174, 175, 165, 176, 175, 176]
a2 = unique(a)

@test SeqIntSet(a1) == SeqIntSet(a2)
