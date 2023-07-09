
# Borrowed from Gridap.jl
macro publish(mod, name)
    quote
      using MMJMesh.$mod: $name; export $name
    end
end


# -------------------------------------------------------------------------------------------------
# Topology
# -------------------------------------------------------------------------------------------------

# ConnectivityList
@publish Topologies ConnectivityList
@publish Topologies push!
@publish Topologies Topology
@publish Topologies nentities
@publish Topologies entity
@publish Topologies nlinks
@publish Topologies links
@publish Topologies link

# Topology
@publish Topologies dim
@publish Topologies addlinks!


