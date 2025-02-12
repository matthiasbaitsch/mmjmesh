module Symbols

using MakieCore
using LinearAlgebra

using MMJMesh
using MMJMesh.MMJBase

export PSymbol, SymbolCollection, symbolfor, stylefor, drawsymbols!, addsymbol!, drawsymbols!


# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

"""
    _bestindex(ads, sds)

Given actual direction vectors `ads` and symbol direction vectors `sds`,
computes the index `i` such that `ads[i]` is the best orientation for the
symbol.
"""
function _bestindex(ads::RealVecVec, sds::RealVecVec)
    badness(a, b) = -a ⋅ b / (norm(a) * norm(b))
    bs = [sum([badness(sd, ad) for ad = ads]) for sd = sds]
    return findmin(bs)[2]
end


struct PSymbol{T}
    marker
    directions::RealVecVec
    angles::RealVec
end

represents(::Type{PSymbol{T}}) where {T} = T
represents(s::PSymbol) = represents(typeof(s))

function symbolfor(o) end
function stylefor(o) end


# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------


struct PlacementList
    points::RealVecVec
    angles::RealVec
    PlacementList() = new(Vector{Vector{Float32}}(undef, 0), Float32[])
end


const SymbolCollection = Dict{PSymbol,PlacementList}


function addsymbol!(sps::SymbolCollection, object, position::RealVec, directions::RealVecVec)
    s = symbolfor(object)
    !haskey(sps, s) && (sps[s] = PlacementList())
    sp = sps[s]
    push!(sp.points, position)
    push!(sp.angles, s.angles[_bestindex(directions, s.directions)])
end


function drawsymbols!(ax, sps::SymbolCollection)
    for (symbol, placements) = sps
        if !isempty(placements.points)
            MakieCore.scatter!(
                ax,
                placements.points |> tomatrix;
                marker=symbol.marker,
                rotation=placements.angles, stylefor(represents(symbol))...
            )
        end
    end
end


# -------------------------------------------------------------------------------------------------
# Symbol collections for special purposes
# -------------------------------------------------------------------------------------------------

# include("symbols.structural2d.jl")

end