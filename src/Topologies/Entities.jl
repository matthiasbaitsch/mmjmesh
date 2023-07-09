function makeentity(d::Int, nn::Int, links::Dict{Tuple{Int, Int}, Vector{Vector{Int}}})
    t = Topology(d, nn)
    for ((from, to), lls) in links
        addlinks!(t, from, to, ConnectivityList(lls))
    end
    return t
end

const L2 = makeentity(1, 2, Dict((1, 0) => [[1, 2]]))

const L3 = makeentity(1, 3, Dict((1, 0) => [[1, 2, 3]]))

const TRI3 = makeentity(2, 3, Dict(
    (2, 0) => [[1, 2, 3]],
    (2, 1) => [[1, 2, 3]],
    (1, 0) => [[1, 2], [2, 3], [1, 3]]
))

const Q4 = makeentity(2, 4, Dict(
    (2, 0) => [[1, 2, 3, 4]],
    (2, 1) => [[1, 2, 3, 4]],
    (1, 0) => [[1, 2], [2, 3], [4, 3], [1, 4]]
))

const TD = Dict(
    (1, 2) => L2,
    (1, 3) => L3,
    (2, 3) => TRI3,
    (2, 4) => Q4
)

entitytopology(d::Int, nn::Int) = TD[(d, nn)]
