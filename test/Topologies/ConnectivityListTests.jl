module ConnectivityListTests

using MMJMesh.Topologies: ConnectivityList, inverse
using Test

# Testee
cl = ConnectivityList([[1, 3, 2], [1, 8, 4, 2], [1, 9, 3, 2], [3, 2]])

# Basic
@test length(cl) == 4
@test length(cl, 2) == 4
@test cl[3] == [1, 9, 3, 2]
@test cl[3, 2] == 9
@test [1, 3, 2] ∈ cl
@test [3, 2] ∈ cl
@test [1, 2] ∉ cl

# Transpose
cltrans = cl'
@test length(cltrans) == 9
@test cltrans[1] == [1, 2, 3]
@test cltrans[2] == [1, 2, 3, 4]
@test cltrans[7] == []
@test cltrans[8] == [2]
@test cltrans[9] == [3]

# Inverse
clinv = inverse(cl)
@test clinv[Set([1, 3, 2])] == 1
@test clinv[Set([1, 9, 3, 2])] == 3
@test clinv[Set([3, 2])] == 4

# Comparison
cl1 = ConnectivityList([[1, 2, 5, 4], [2, 3, 5], [3, 6, 5]])
cl2 = ConnectivityList([[1, 2, 5, 4], [2, 3, 5], [3, 6, 5]])
@test cl1 == cl2


end # module