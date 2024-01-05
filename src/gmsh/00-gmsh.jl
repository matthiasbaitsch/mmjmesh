module Gmsh

# Modules needed by this module
using Lerche
using Tables
using PrettyTables
import MMJMesh.Utilities: SeqIntSet

# Parts
include("objecttable.jl")
include("arrayscanner.jl")
include("gmshmesh.jl")
include("readmesh.jl")

end