"""
 More or less verbatim translation of msh file format to Julia structs.
"""

struct MeshFormat
    version::Float64
    fileType::Int
    dataSize::Int
end


struct PhysicalName
    dimension::Int
    tag::Int
    name::String
end

Base.show(io::IO, mf::MeshFormat) = println(io, "# MeshFormat\n\n  version: $(mf.version)\n filetype: $(mf.fileType)\n datasize: $(mf.dataSize)")
Base.show(io::IO, pn::PhysicalName) = println(io, "PhysicalName[dimension=$(pn.dimension), tag=$(pn.tag), name=$(pn.name)]")

struct PhysicalNameCollection
    nNames::Int
    names::Vector{PhysicalName}
end

const PT_CONF = set_pt_conf(crop=:horizontal)

function Base.show(io::IO, pnc::PhysicalNameCollection)
    if length(pnc.names) == 0
        println("# No physical names")
    else
        println("# Physical names\n")
        pretty_table_with_conf(PT_CONF, io, ObjectTable(pnc.names))
    end
end

Base.length(ebc::PhysicalNameCollection) = ebc.nNames
Base.getindex(ebc::PhysicalNameCollection, idx) = ebc.names[idx]

struct BoundingBox
    bounds::Vector{Float64}
end

struct Point
    tag::Int
    position::Vector{Float64}
    physicalTags::Vector{Int}
end

struct Entity
    tag::Int
    boundingBox::BoundingBox
    physicalTags::Vector{Int}
    boundingEntities::Vector{Int}
end



struct EntityCollection
    # Note that entity numbering does not
    # always start at 1, see e.g. complex-g1.msh
    points::Dict{Int,Any}
    curves::Dict{Int,Any}
    surfaces::Dict{Int,Any}
    volumes::Dict{Int,Any}
    allEntities::Dict{Int,Any}
    EntityCollection(points, curves, surfaces, volumes) = begin
        toDict(points) = reduce((d, e) -> (d[e.tag] = e; d), points, init=Dict{Int,Any}())
        pd = toDict(points)
        cd = toDict(curves)
        sd = toDict(surfaces)
        vd = toDict(volumes)
        return new(
            pd, cd, sd, vd,
            Dict([0 => pd, 1 => cd, 2 => sd, 3 => vd])
        )
    end
end
Base.getindex(ec::EntityCollection, dim::Int) = ec.allEntities[dim]

function Base.show(io::IO, ec::EntityCollection)
    println(io, "# Entities")
    if !isempty(ec.points)
        println(io, "\n## Points\n")
        pretty_table_with_conf(PT_CONF, io, ObjectTable(ec.points))
    end
    if !isempty(ec.curves)
        println(io, "\n## Curves\n")
        pretty_table_with_conf(PT_CONF, io, ObjectTable(ec.curves))
    end
    if !isempty(ec.surfaces)
        println(io, "\n## Surfaces\n")
        pretty_table_with_conf(PT_CONF, io, ObjectTable(ec.surfaces))
    end
    if !isempty(ec.volumes)
        println(io, "\n## Volumes\n")
        pretty_table_with_conf(PT_CONF, io, ObjectTable(ec.volumes))
    end
end

abstract type Block end
abstract type BlockCollection end

struct NodeBlock <: Block
    entityDim::Int
    entityTag::Int
    parametric::Bool
    nodeTags::Vector{Int}
    coordinates::Matrix{Float64}
end

struct NodeBlockCollection <: BlockCollection
    nBlocks::Int
    nNodes::Int
    minNodeTag::Int
    maxNodeTag::Int
    blocks::Vector{NodeBlock}
end

Base.length(ebc::NodeBlockCollection) = ebc.nBlocks
Base.getindex(ebc::NodeBlockCollection, i) = ebc.blocks[i]

function Base.show(io::IO, ec::NodeBlockCollection)
    println(io, "# Node Blocks\n")
    println(io, "    nBlocks: ", ec.nBlocks)
    println(io, "     nNodes: ", ec.nNodes)
    println(io, " minNodeTag: ", ec.minNodeTag)
    println(io, " maxNodeTag: ", ec.maxNodeTag)
    if !isempty(ec.blocks)
        println(io)
        pretty_table_with_conf(PT_CONF, io, ObjectTable(ec.blocks))
    end
end

struct ElementBlock <: Block
    entityDim::Int
    entityTag::Int
    elementType::Int
    elementTags::Vector{Int}
    nodeTags::Matrix{Int}
end

struct ElementBlockCollection <: BlockCollection
    nBlocks::Int
    nElements::Int
    minElementTag::Int
    maxElementTag::Int
    blocks::Vector{ElementBlock}
end

function Base.show(io::IO, ec::ElementBlockCollection)
    println(io, "# Element Blocks\n")
    println(io, "       nBlocks: ", ec.nBlocks)
    println(io, "     nElements: ", ec.nElements)
    println(io, " minElementTag: ", ec.minElementTag)
    println(io, " maxElementTag: ", ec.maxElementTag)
    if !isempty(ec.blocks)
        println(io)
        pretty_table_with_conf(PT_CONF, io, ObjectTable(ec.blocks))
    end
end

Base.length(ebc::ElementBlockCollection) = ebc.nBlocks
Base.getindex(ebc::ElementBlockCollection, i) = ebc.blocks[i]
blocksByDimension(bc::BlockCollection, dim::Int) = filter(b -> b.entityDim == dim, bc.blocks)

function nodeTags(ebc::ElementBlockCollection, dim::Int)
    nts = []
    for eb ∈ ebc.blocks
        if eb.entityDim == dim
            if length(nts) == 0
                nts = eb.nodeTags
            else
                nts = vcat(nts, eb.nodeTags)
            end
        end
    end
    return nts
end

struct GmshMesh
    meshFormat::MeshFormat
    physicalNames::PhysicalNameCollection
    entities::EntityCollection
    nodeBlocks::NodeBlockCollection
    elementBlocks::ElementBlockCollection
end

function Base.show(io::IO, m::GmshMesh)
    println(io, m.meshFormat)
    println(io, m.physicalNames)
    println(io, m.entities)
    println(io, m.nodeBlocks)
    println(io, m.elementBlocks)
end

dimension(m::GmshMesh) = maximum(n -> n.entityDim, m.elementBlocks.blocks)

function blockName(m::GmshMesh, b::Block)
    entity = m.entities[b.entityDim][b.entityTag]
    physicalTags = entity.physicalTags
    if length(physicalTags) == 0
        return "_$(b.entityDim).$(b.entityTag)"
    elseif length(physicalTags) == 1
        return m.physicalNames[physicalTags[1]].name
    else
        error("Should not happen")
    end
end

function getnodes(m::GmshMesh)
    Nn = m.nodeBlocks.nNodes
    nodes = zeros(Nn, 3)
    for nb ∈ m.nodeBlocks.blocks
        nodes[nb.nodeTags, :] = nb.coordinates
    end
    if norm(nodes[:, 3]) == 0
        nodes = nodes[:, 1:2]
    end
    return nodes'
end

