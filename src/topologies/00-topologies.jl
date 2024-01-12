module Topologies

# Modules needed by this module
using Printf
using MMJMesh

# Exports
## ConnectivityList
export ConnectivityList
export inverse, maxlinkssize
## Topology
export Topology
export addlinks!, nentities, entity, dimension, isanonymous, nlinks, links

# Parts
include("connectivitylist.jl")
include("topology.jl")
include("entitytopologies.jl")

end