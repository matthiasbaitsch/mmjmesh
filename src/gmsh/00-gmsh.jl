module Gmsh

# Modules needed by this module
using Lerche
using Tables
using LinearAlgebra
using PrettyTables

using MMJMesh.Topologies
import MMJMesh.Meshes: Mesh, populatepredfinedgroups!
import MMJMesh.MMJBase: SeqIntSet
import MMJMesh.Groups: EntityGroup

# Exports
export Mesh

# Parts
include("objecttable.jl")
include("arrayscanner.jl")
include("gmshmesh.jl")
include("readmesh.jl")

# Create Mesh
function Mesh(filepath::String)
    gm = readmesh(filepath)
    D = dimension(gm)
    m = Mesh(coordinates(gm), D)

    # Collect empty groups
    groupnamesbytag = Dict{Int,Symbol}()
    for g ∈ gm.physicalnames.names
        name = Symbol(g.name)
        m.groups[name] = EntityGroup(g.dimension, Int[])
        groupnamesbytag[g.tag] = name
    end

    # Collect entities
    for eb ∈ gm.elementblocks.blocks
        d0 = eb.entitydim

        if d0 > 0
            # Links
            indexes = addlinks!(m.topology, d0, 0, eb.tags, ConnectivityList(eb.nodetags))

            # Groups
            for tag ∈ gm.entities[d0][eb.entityTag].physicaltags
                name = groupnamesbytag[tag]
                m.groups[name] = m.groups[name] ∪ EntityGroup(d0, indexes)
            end
        end
    end

    # Make sure that all links are created
    for d1 ∈ 1:D-1
        links(m.topology, D, d1)
    end

    # Predefined groups, make sure that all links are created
    populatepredfinedgroups!(m)

    return m
end

end