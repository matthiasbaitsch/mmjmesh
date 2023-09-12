
# Borrowed from Gridap.jl
macro publish(mod, name)
  quote
    using MMJMesh.$mod: $name
    export $name
  end
end


# -------------------------------------------------------------------------------------------------
# Topologies
# -------------------------------------------------------------------------------------------------

# ConnectivityList
# @publish Topologies ConnectivityList
# @publish Topologies entity
# @publish Topologies link
# @publish Topologies links
# @publish Topologies nentities
# @publish Topologies nconnections
# @publish Topologies maxlinksize
# @publish Topologies push!
# @publish Topologies Topology

# # Topology
# @publish Topologies addlinks!
# @publish Topologies dimension


# -------------------------------------------------------------------------------------------------
# Geometries
# -------------------------------------------------------------------------------------------------

# Point
# @publish Geometries boundingbox
# @publish Geometries center
# @publish Geometries coordinates
# @publish Geometries pdim
# @publish Geometries gdim
# @publish Geometries measure
# @publish Geometries Point

# # Geometry
# @publish Geometries coordinates
# @publish Geometries Geometry
# @publish Geometries squeeze!

# Utilities
# @publish Utilities makemeshonrectangle
# @publish Utilities plotMesh
# @publish Utilities PlotMeshConfiguration


# -------------------------------------------------------------------------------------------------
# Meshes
# -------------------------------------------------------------------------------------------------

# @publish Meshes Mesh

