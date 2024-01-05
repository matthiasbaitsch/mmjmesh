using Test

import MMJMesh.Utilities: SeqIntSet

s = SeqIntSet([1, 10, -50, 2, 11, 9, 8, 0, 100])
@test length(s.sequences) == 4
@test length(s) == 9
@test collect(s) == [-50, 0, 1, 2, 8, 9, 10, 11, 100]
@test -1000 ∉ s
@test -51 ∉ s
@test 5 ∉ s
@test 12 ∉ s
@test 101 ∉ s
@test 1000 ∉ s
@test all([i ∈ s for i ∈ s])

s = SeqIntSet([1])
@test length(s) == 1
@test collect(s) == [1]
@test 1 ∈ s
@test -1 ∉ s
@test 2 ∉ s

s = SeqIntSet(Vector{Int}())
@test isempty(s)
@test collect(s) == Vector{Int}()
@test 5 ∉ s
