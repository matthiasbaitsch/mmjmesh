module Gmsh

# Modules needed by this module
using Gmsh
using Lerche
using Tables
using LinearAlgebra
using PrettyTables

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.MMJBase
using MMJMesh.Topologies


# Exports

# Parts
include("objecttable.jl")
include("arrayscanner.jl")
include("gmshmesh.jl")
include("readmesh.jl")

function _readmsh(filepath)
    gm = readmesh(filepath)
    D = dimension(gm)
    m = Mesh(coordinates(gm), D)

    ng(name) = Symbol(string(name) * "0")

    # Collect empty groups
    groupnamesbytag = Dict{Int,Symbol}()
    for g ∈ gm.physicalnames.names
        name = Symbol(g.name)
        definegroup!(m, g.dimension, name, Int[])
        groupnamesbytag[g.tag] = name
        if g.dimension == 1
            definegroup!(m, 0, ng(name), Int[])
        end
    end

    # Collect entities
    for eb ∈ gm.elementblocks.blocks
        d0 = eb.entitydim

        if d0 > 0
            # Links
            indices = addlinks!(m.topology, d0, 0, eb.tags, ConnectivityList(eb.nodetags))

            # Groups
            for tag ∈ gm.entities[d0][eb.entityTag].physicaltags
                name = groupnamesbytag[tag]
                m.groups[name] = m.groups[name] ∪ Group{MeshEntity{d0}}(indices)

                if d0 == 1
                    nn = ng(name)
                    m.groups[nn] = m.groups[nn] ∪ Group{MeshEntity{0}}(reshape(eb.nodetags, :))
                end
            end
        end
    end

    # Make sure that all links are created
    for d1 ∈ 1:D-1
        links(m.topology, D, d1)
    end

    return m
end

function _readgeo(filepath; verbosity)
    try
        f = tempname() * ".msh"
        gmsh.initialize()
        gmsh.option.setNumber("General.Verbosity", verbosity)
        gmsh.open(filepath)
        gmsh.write(f)
        gmsh.finalize()
        m = _readmsh(f)
        rm(f)
        return m
    catch e
        gmsh.finalize()
        throw(e)
    end
end

function _readgeostring(s; verbosity)
    f = tempname() * ".geo"
    io = open(f, "a")
    println(io, s)
    close(io)
    m = _readgeo(f, verbosity=verbosity)
    rm(f)
    return m
end

# Create Mesh
function MMJMesh.Meshes.Mesh(filepath::String; verbosity=0)
    if endswith(filepath, ".msh")
        return _readmsh(filepath)
    end
    if endswith(filepath, ".geo")
        return _readgeo(filepath; verbosity=verbosity)
    end
    return _readgeostring(filepath; verbosity=verbosity)
end

end