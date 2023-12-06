# -------------------------------------------------------------------------------------------------
# Interval
# -------------------------------------------------------------------------------------------------
@test Interval(1, 2).a == 1
@test Interval(1, 2).b == 2
@test dim(Interval(1, 2)) == 1
@test in(Interval(1, 2), 1)
@test in(Interval(1, 2), 1.5)
@test in(Interval(1, 2), 2)
@test !in(Interval(1, 2), 0.9)
@test !in(Interval(1, 2), 2.1)
@test in(Interval(1, 2), 0.9, eps=0.1)
@test in(Interval(1, 2), 2.1, eps=0.1)
@test discretize(Interval(1, 2), 3) == [1.0, 1.5, 2.0]
