"""
    Topology(entities, links)

A `Topology` stores links between entities of a mesh.

Entities have a (parametric/topological) dimension. We use the terminology
- d=0: Node
- d=1: Edge
- d=2: Face
- d=3: Region
although these terms do not show in the code very much.

When constructing a topology, entities are added by specifying their IDs and 
links to vertexes. Each entity of dimension d > 1 define additional entities
of dimensions d-1, ..., 1. If not specified explicitly, these entities are
added as needed but no not have IDs. Thus, they are called anonymous. 

Remarks:
- Links refer to the index of an entity not the ID.
- IDs are mostly for compatibility with software which relies on numbering.
"""
struct Topology{D}
    entities::Dict{Int,Vector{Int}}
    links::Dict{Tuple{Int,Int},ConnectivityList}
end

function Topology(d::Int, nn::Int)
    return Topology{d}(Dict(0 => collect(1:nn)), Dict{Tuple{Int,Int},ConnectivityList}())
end

dim(t::Topology{D}) where {D} = D
nentities(t::Topology, d::Int) = haskey(t.entities, d) ? length(t.entities[d]) : 0
entity(t::Topology, dim::Int, idx::Int) = t.entities[dim][idx]

"""
    addlinks!(t, d0, d1, [ids,] cl)

Add links to the topology.
"""
function addlinks!(t::Topology, d0::Int, d1::Int, ids::Vector{Int}, cl::ConnectivityList)
    t.entities[d0] = ids
    t.links[(d0, d1)] = cl
    return nothing
end

addlinks!(t::Topology, d0::Int, d1::Int, cl::ConnectivityList) = addlinks!(t, d0, d1, collect(1:length(cl)), cl)
addlinks!(t::Topology, d0::Int, d1::Int, links::Vector{Vector{Int}}) = addlinks!(t, d0, d1, ConnectivityList(links))
nlinks(t::Topology, d0::Int, d1::Int) = haskey(t.links, (d0, d1)) ? length(t.links[(d0, d1)]) : 0
link(t::Topology, d0::Int, d1::Int, idx::Int) = t.links[(d0, d1)][idx]

function links(t::Topology{D}, d0::Int, d1::Int) where {D}
    key = (d0, d1)

    if (haskey(t.links, key))
        return t.links[key]
    elseif d0 > 0 && d1 == 0                                        # a.1) Node links
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
    elseif d0 > d1 && d1 > 0                                        # a.2) Entity links
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
    elseif d0 == d1                                                 #a.3) Same-same links
        ne = d0 == 0 ? nentities(t, 0) : length(links(t, d0, 0))
        linkslist = [Int[] for _ in 1:ne]
        if d0 > 0
            for ll ∈ links(t, d0 - 1, d0), i ∈ ll, j ∈ ll
                if i != j
                    push!(linkslist[i], j)
                end
            end
        end
        cl = ConnectivityList(linkslist)        
    elseif d1 > d0                                                  # a.4) Upward links
        cl = links(t, d1, d0)'
    end

    # Cache
    t.links[key] = cl

    # Return
    return cl
end

function Base.show(io::IO, t::Topology; all::Bool=false)
    println(io, "*-----------------------------------------------------------------------------------------------------------*")
    println(io, "Topology{", dim(t), "}")
    println(io, "-------------------------------------------------------------------------------------------------------------")
    println(io, "Entities")
    for d in 0:dim(t)
        if haskey(t.entities, d)
            println(io, "$d: $(t.entities[d])")
        end
    end

    println(io, length(t.links) == 0 ? "No links" : "Links")
    for ((d0, d1), cl) in t.links    
        if all || (d0 ≠ d1 && haskey(t.entities, d0) && haskey(t.entities, d1))
            print(io, "$d0 -> $d1\n$cl")
        end
    end
    println(io, "*-----------------------------------------------------------------------------------------------------------*\n")
end

Base.show(t::Topology; all::Bool=false) = show(stdout, t, all=all)