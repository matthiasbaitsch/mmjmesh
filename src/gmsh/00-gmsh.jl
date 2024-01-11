module Gmsh

# Modules needed by this module
using Lerche
using Tables
using LinearAlgebra
using PrettyTables
import MMJMesh.Meshes: Mesh
import MMJMesh.MMJBase: SeqIntSet

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
    return Mesh(coordinates(gm), nodetags(gm.elementblocks, D), D)
end

end