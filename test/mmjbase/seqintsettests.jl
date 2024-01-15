using Test

using MMJMesh
import MMJMesh.MMJBase: SeqIntSet, nsequences, start, step, nsteps, stop


# -------------------------------------------------------------------------------------------------
# Constructor
# -------------------------------------------------------------------------------------------------

# Three sequences
a = [1, 3, 5, 6, 7, 8, 9, 11, 14, 17]
set = SeqIntSet(a);
@test set.index == [1, 4, 8, 11]
@test set.start == [1, 6, 11]
@test set.step == [2, 1, 3]

# No sequences
a = [1, 2, 4, 5, 7]
set = SeqIntSet(a);
@test set.start == a
@test set.step == [1, 2, 1, 2, 0]
@test set.index == collect(1:(length(a)+1))

# One sequence
a = [3, 5, 7, 9, 11]
set = SeqIntSet(a);
@test set.start == [3]
@test set.step == [2]
@test set.index == [1, 6]

# One sequence, one single number at the end
a = [3, 5, 7, 9, 11, 100]
set = SeqIntSet(a);
@test set.start == [3, 100]
@test set.step == [2, 0]
@test set.index == [1, 6, 7]

# Two sequences, three single numbers
a = [0, 2, 3, 4, 6, 8, 10, 30, 31]
set = SeqIntSet(a);
@test set.index == [1, 2, 5, 8, 9, 10]
@test set.start == [0, 2, 6, 30, 31]
@test set.step == [2, 1, 2, 1, 0]


# -------------------------------------------------------------------------------------------------
# Methods
# -------------------------------------------------------------------------------------------------

function validate(a::AbstractVector, set::SeqIntSet)
    b = unique(sort(a))
    for i in a
        @test i ∈ set
    end
    for i in set
        @test i ∈ a
        @test i ∈ set
    end
    for i ∈ setdiff(first(set):last(set), a)
        @test i ∉ set
    end
    for i ∈ eachindex(b)
        @test set[i] == b[i]
    end
    @test set == b
    @test collect(set) == b
end


# Methods basic tests
a = [0, 2, 3, 4, 6, 8, 10, 30, 31]
set = SeqIntSet(a)

@test nsequences(set) == 5
@test first(set) == 0
@test length(set) == 9
@test size(set) == (9,)
@test !isempty(set)
@test last(set) == 31
@test start(set, 1) == 0
@test step(set, 1) == 2
@test nsteps(set, 1) == 0
@test stop(set, 1) == 0


# Empty
set = SeqIntSet(Vector{Int}());
@test isempty(set)
@test collect(set) == Vector{Int}()
@test 5 ∉ set


# One element
set = SeqIntSet([1])
@test length(set) == 1
@test collect(set) == [1]
@test 1 ∈ set
@test -1 ∉ set
@test 2 ∉ set
validate([1], set)


# Many elements
a = [1, 10, -50, 2, 11, 9, 8, 0, 100]
set = SeqIntSet(a);

@test nsequences(set) == 4
@test size(set) == (9,)
@test length(set) == 9
@test collect(set) == [-50, 0, 1, 2, 8, 9, 10, 11, 100]
@test -1000 ∉ set
@test -51 ∉ set
@test 1 ∈ set
@test 5 ∉ set
@test 12 ∉ set
@test 100 ∈ set
@test 101 ∉ set
@test 1000 ∉ set
@test set[1] == -50
@test set[2] == 0
@test set[3] == 1
@test set[4] == 2
@test set[5] == 8
@test set[8] == 11
@test set[9] == 100
validate(a, set)


# Test case from an actual mesh
a1 = [1, 2, 1, 12, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 22, 12, 23, 22, 33, 23, 34, 33, 44, 34, 45, 44, 55, 45, 56, 55, 66, 56, 67, 66, 77, 67, 78, 77, 88, 78, 89, 88, 99, 89, 100, 99, 110, 100, 111, 110, 121, 111, 122, 121, 132, 122, 133, 132, 143, 133, 144, 143, 154, 144, 155, 154, 165, 166, 167, 155, 166, 167, 168, 168, 169, 169, 170, 170, 171, 171, 172, 172, 173, 173, 174, 174, 175, 165, 176, 175, 176]

validate(a1, SeqIntSet(a1))
@test SeqIntSet(a1) == SeqIntSet(unique(a1))


# Set operations
s1 = SeqIntSet([1, 2, 3, 5, 9, 10])
s2 = SeqIntSet([1, 2, 3, 4, 5, 6])

@test s1 ∪ s2 == [1, 2, 3, 4, 5, 6, 9, 10]
@test s1 ∩ s2 == [1, 2, 3, 5]
@test setdiff(s1, s2) == [9, 10]