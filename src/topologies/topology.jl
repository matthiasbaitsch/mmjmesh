"""
    Topology(entities, links)

A `Topology` stores connections between entities of a mesh. The approach
is based on the mesh implementation in the FEniCSx project.

Mesh entities have a parametric dimension ``d``. We use the terminology

- ``d=0``: Node
- ``d=1``: Edge
- ``d=2``: Face
- ``d=3``: Solid

Connections between entities (of two not necessarily different dimensions) 
are stored in `ConnectivityList`s. A `ConnectivityList` contains links
which are arrays of indexes of other mesh entities. Hence, in the code, 
a connectivity list is called `linkslist`, one element in the list is called `links`.

When constructing a topology, entities are added by specifying their IDs and 
connectivities to vertexes. Each entity of dimension d > 1 defines additional
entities of dimensions d-1, ..., 1 (such as a face comes with edges). If not 
specified explicitly, these entities are added as needed but do not have IDs. 
Thus, they are called anonymous. 

Remarks:
- Links refer to the index of an entity not the ID.
- IDs are mostly for compatibility with software which relies on numbering.

Example (anonymous entities are in parentheses):

        4   (3)   5   (7)   6
        *---------*---------*
        |         |         |
     (4)|    1    |(2) 2    |(6)
        |         |         |
        *---------*---------*
        1   (1)   2   (5)   3

- Entities
    - 0: [1, 2, 3, 4, 5, 6]
    - 2: [1, 2]
- Links
    - (2, 0): [[1, 2, 5, 4], [2, 3, 6, 5]]
- Generated and upward links
    - (2, 2): [[2], [1]]
    - (2, 1): [[1, 2, 3, 4], [5, 6, 7, 2]]
    - (1, 2): [[1], [1, 2], [1], [1], [2], [2], [2]]
    - (1, 1): [[2, 4, 5], [1, 3, 5, 7], [2, 4, 7], [1, 3], [1, 2, 6], [5, 7], [2, 3, 6]]
    - (1, 0): [[1, 2], [2, 5], [4, 5], [1, 4], [2, 3], [3, 6], [5, 6]]
    - (0, 2): [[1], [1, 2], [2], [1], [1, 2], [2]]
    - (0, 1): [[1, 4], [1, 2, 5], [5, 6], [3, 4], [2, 3, 7], [6, 7]]
    - (0, 0): [[], [], [], [], [], []]
"""
struct Topology{D}
    entities::Dict{Int,Vector{Int}}
    links::Dict{Tuple{Int,Int},ConnectivityList}
end


"""
    Topology(d, nn)

Construct topology of parametric dimension `d` with `nn` nodes.
"""
function Topology(d::Int, nn::Int)
    return Topology{d}(Dict(0 => collect(1:nn)), Dict{Tuple{Int,Int},ConnectivityList}())
end

dimension(::Topology{D}) where {D} = D
isanonymous(t::Topology, d::Int) = !haskey(t.entities, d)

# -------------------------------------------------------------------------------------------------
# Entities
# -------------------------------------------------------------------------------------------------

function nentities(t::Topology{D}, d::Int, anonymous::Bool=true) where {D}
    if d <= D
        if haskey(t.entities, d)
            return length(t.entities[d])
        elseif anonymous
            return length(links(t, d, 0))
        else
            return 0
        end
    else
        return 0
    end
end

entity(t::Topology, dim::Int, idx::Int) = t.entities[dim][idx]

# -------------------------------------------------------------------------------------------------
# Links
# -------------------------------------------------------------------------------------------------

"""
    addlinks!(t::Topology, d0::Int, d1::Int, ids::Vector{Int}, cl::ConnectivityList)

Add links from entities of dimension `d0` to entities of dimension `d1` by specifying the
entities `ids` and the `ConnectivityList` `cl`.
"""
function addlinks!(t::Topology{D}, d0::Int, d1::Int, ids::Vector{Int}, cl::ConnectivityList) where {D}
    @assert d0 <= D && d0 > d1 "Only downward links with d0 < $D are allowed."
    t.entities[d0] = ids
    t.links[(d0, d1)] = cl
    return nothing
end

"""
    addlinks!(t::Topology, d0::Int, d1::Int, cl)

Add links from entities of dimension `d0` to entities of dimension `d1` by specifying the connectivity list.
Entity IDs are generated automatically. The parameter `cl` can be a `ConnectivityList` or a vector of integer vectors.
"""
addlinks!(t::Topology, d0::Int, d1::Int, cl::ConnectivityList) = addlinks!(t, d0, d1, collect(1:length(cl)), cl)
addlinks!(t::Topology, d0::Int, d1::Int, cl::Vector{Vector{Int}}) = addlinks!(t, d0, d1, ConnectivityList(cl))

"""
    links(t::Topology, d0::Int, d1::Int)

Links from entities of dimension `d0` to entities of dimension `d1` in the form of
a `ConnectivityList`.
"""
function links(t::Topology{D}, d0::Int, d1::Int) where {D}
    @assert d0 <= D && d1 <= D
    key = (d0, d1)

    if (haskey(t.links, key))
        return t.links[key]
    elseif d0 > 0 && d1 == 0                                        # a.1) Links to nodes
        visited = Set{Set{Int}}()
        linkslist = Vector{Vector{Int}}()
        for gl ∈ links(t, D, 0)
            for ll ∈ links(entitytopology(D, length(gl)), d0, 0)
                idxs = gl[ll]
                idxsset = Set(idxs)
                if idxsset ∉ visited
                    push!(visited, idxsset)
                    push!(linkslist, idxs)
                end
            end
        end
        cl = ConnectivityList(linkslist)
    elseif d0 > d1 && d1 > 0                                        # a.2) Links entities of lower dimension
        im = inverse(links(t, d1, 0))
        linkslist = Vector{Vector{Int}}()
        for gl ∈ links(t, d0, 0)
            llinks = Vector{Int}()
            for idxs ∈ links(entitytopology(D, length(gl)), d1, 0)
                push!(llinks, im[Set(gl[idxs])])
            end
            push!(linkslist, llinks)
        end
        cl = ConnectivityList(linkslist)
    elseif d0 == d1                                                 #a.3) Links to entities of same dimension
        ne = d0 == 0 ? nentities(t, 0) : length(links(t, d0, 0))
        linkslist = [Int[] for _ in 1:ne]
        if d0 > 0
            for ll ∈ links(t, d0 - 1, d0), i ∈ ll, j ∈ ll
                if i != j
                    push!(linkslist[i], j)
                end
            end
        end
        foreach(sort!, linkslist)
        cl = ConnectivityList(linkslist)
    elseif d1 > d0                                                  # a.4) Upward links
        cl = links(t, d1, d0)'
    end

    t.links[key] = cl
    return cl
end


# -------------------------------------------------------------------------------------------------
# IO
# -------------------------------------------------------------------------------------------------


"""
    show([io, ]t::Topology[, all=false])

Show the topology `t` on `io` (default: `stdout`). If `all` is `true`, anonymous links are
shown as well.
"""
function Base.show(io::IO, t::Topology{D}; all::Bool=false) where {D}
    # XXX does not work 
    cio = IOContext(io, :short => true)
    println(io, "Topology{$D}")
    println(io, "Entities")
    for d in 0:dimension(t)
        if haskey(t.entities, d)
            println(cio, "         $d: $(t.entities[d])")
        end
    end
    if length(t.links) > 0
        processed = Set()
        println(io, "Links")
        for d0 ∈ D:-1:0, d1 ∈ d0-1:-1:0
            if !isanonymous(t, d0) && !isanonymous(t, d1)
                println(cio, "    ($d0, $d1): $(links(t, d0, d1))")
                push!(processed, (d0, d1))
            end
        end
        if all
            println("Generated and upward links")
            for d0 ∈ D:-1:0, d1 ∈ D:-1:0
                if (d0, d1) ∉ processed
                    println(cio, "    ($d0, $d1): $(links(t, d0, d1))")
                end
            end
        end
    else
        println("No links")
    end
end

Base.show(t::Topology; all::Bool=false) = show(stdout, t, all=all)