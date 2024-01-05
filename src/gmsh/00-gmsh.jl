module Gmsh

# Modules needed by this module
using Lerche
using Tables
using PrettyTables

# Parts
include("objecttable.jl")
include("arrayscanner.jl")
include("gmshmesh.jl")
include("readmesh.jl")

end