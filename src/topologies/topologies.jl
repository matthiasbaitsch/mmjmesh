module Topologies


# Topology
export Topology
export addlinks!
export entity
export isanonymous
export links
export nentities

# ConnectivityList
export ConnectivityList
export inverse
export maxlinksize

using Printf


include("connectivitylist.jl")

include("topology.jl")

include("entitytopologies.jl")

end