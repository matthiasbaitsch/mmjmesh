module Topologies

# Exports
## ConnectivityList
export ConnectivityList
export inverse
export maxlinkssize
## Topology
export Topology
export addlinks!
export entity
export isanonymous
export links
export nentities

# Modules needed by this module
using Printf

# Parts
include("connectivitylist.jl")
include("topology.jl")
include("entitytopologies.jl")

end