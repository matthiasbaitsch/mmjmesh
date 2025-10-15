# More or less verbatim translation of Gmsh's .msh file format to Julia structs

struct MeshFormat
    version::Float64
    filetype::Int
    datasize::Int
end

struct PhysicalName
    dimension::Int
    tag::Int
    name::String
end

struct PhysicalNameCollection
    nnames::Int
    names::Vector{PhysicalName}
end

struct BoundingBox
    bounds::Vector{Float64}
end

struct Point
    tag::Int
    position::Vector{Float64}
    physicaltags::Vector{Int}
end

struct Entity
    tag::Int
    boundingbox::BoundingBox
    physicaltags::Vector{Int}
    boundingentities::Vector{Int}
end

struct EntityCollection
    # Note that entity numbering does not
    # always start at 1, see e.g. complex-g1.msh
    points::Dict{Int,Any}
    curves::Dict{Int,Any}
    surfaces::Dict{Int,Any}
    volumes::Dict{Int,Any}
    allentities::Dict{Int,Any}
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

abstract type Block end

abstract type BlockCollection end

struct NodeBlock <: Block
    entitydim::Int
    entityTag::Int
    parametric::Bool
    tags::SeqIntSet
    coordinates::Matrix{Float64}
end

struct NodeBlockCollection <: BlockCollection
    nblocks::Int
    nnodes::Int
    minnodetag::Int
    maxnodetag::Int
    blocks::Vector{NodeBlock}
end

struct ElementBlock <: Block
    entitydim::Int
    entityTag::Int
    type::Int
    tags::SeqIntSet
    nodetags::Matrix{Int}
end

struct ElementBlockCollection <: BlockCollection
    nblocks::Int
    nelements::Int
    minelementtag::Int
    maxelementtag::Int
    blocks::Vector{ElementBlock}
end

struct GmshMesh
    meshFormat::MeshFormat
    physicalnames::PhysicalNameCollection
    entities::EntityCollection
    nodeBlocks::NodeBlockCollection
    elementblocks::ElementBlockCollection
end


# MeshFormat and PhysicalName methods

Base.show(io::IO, mf::MeshFormat) = print(io, "# MeshFormat\n\n  version: $(mf.version)\n filetype: $(mf.filetype)\n datasize: $(mf.datasize)\n")
Base.show(io::IO, pn::PhysicalName) = print(io, "PhysicalName[dimension=$(pn.dimension), tag=$(pn.tag), name=$(pn.name)]")
Base.show(io::IO, bb::BoundingBox) = print(io, bb.bounds)


# PhysicalNameCollection methods

function Base.show(io::IO, pnc::PhysicalNameCollection)
    if length(pnc.names) == 0
        println("# No physical names")
    else
        println("# Physical names\n")
        pretty_table(io, ObjectTable(pnc.names))
    end
end

Base.length(ebc::PhysicalNameCollection) = ebc.nnames
Base.getindex(ebc::PhysicalNameCollection, idx) = ebc.names[idx]


# EntityCollection methods

Base.getindex(ec::EntityCollection, dim::Int) = ec.allentities[dim]

function Base.show(io::IO, ec::EntityCollection)
    println(io, "# Entities")
    if !isempty(ec.points)
        println(io, "\n## Points\n")
        pretty_table(io, ObjectTable(ec.points))
    end
    if !isempty(ec.curves)
        println(io, "\n## Curves\n")
        pretty_table(io, ObjectTable(ec.curves))
    end
    if !isempty(ec.surfaces)
        println(io, "\n## Surfaces\n")
        pretty_table(io, ObjectTable(ec.surfaces))
    end
    if !isempty(ec.volumes)
        println(io, "\n## Volumes\n")
        pretty_table(io, ObjectTable(ec.volumes))
    end
end


# NodeBlockCollection methods

Base.length(ebc::NodeBlockCollection) = ebc.nblocks
Base.getindex(ebc::NodeBlockCollection, i) = ebc.blocks[i]

function Base.show(io::IO, ec::NodeBlockCollection)
    println(io, "# Node Blocks\n")
    println(io, "    nblocks: ", ec.nblocks)
    println(io, "     nnodes: ", ec.nnodes)
    println(io, " minnodetag: ", ec.minnodetag)
    println(io, " maxnodetag: ", ec.maxnodetag)
    if !isempty(ec.blocks)
        println(io)
        pretty_table(io, ObjectTable(ec.blocks))
    end
end


# ElementBlockCollection methods

function Base.show(io::IO, ec::ElementBlockCollection)
    println(io, "# Element Blocks\n")
    println(io, "       nblocks: ", ec.nblocks)
    println(io, "     nelements: ", ec.nelements)
    println(io, " minelementtag: ", ec.minelementtag)
    println(io, " maxelementtag: ", ec.maxelementtag)
    if !isempty(ec.blocks)
        println(io)
        pretty_table(io, ObjectTable(ec.blocks))
    end
end

Base.length(ebc::ElementBlockCollection) = ebc.nblocks

function nodetags(ebc::ElementBlockCollection, dim::Int)
    nts = Vector{Vector{Int}}()
    for eb ∈ ebc.blocks        
        if eb.entitydim == dim            
            for i ∈ 1:size(eb.nodetags, 1)
                push!(nts, eb.nodetags[i, :])
            end
        end
    end
    return nts
end


# GmshMesh methods

MMJMesh.dimension(m::GmshMesh) = maximum(n -> n.entitydim, m.elementblocks.blocks)

function blockname(m::GmshMesh, b::Block)
    entity = m.entities[b.entitydim][b.entityTag]
    physicaltags = entity.physicaltags
    if length(physicaltags) == 0
        return "_$(b.entitydim).$(b.entityTag)"
    elseif length(physicaltags) == 1
        return m.physicalnames[physicaltags[1]].name
    else
        error("Should not happen")
    end
end

function coordinates(m::GmshMesh)
    Nn = m.nodeBlocks.nnodes
    nodes = zeros(3, Nn)
    for nb ∈ m.nodeBlocks.blocks
        nodes[:, nb.tags] = nb.coordinates'
    end
    if norm(nodes[3, :]) == 0
        nodes = nodes[1:2, :]
    end
    return nodes
end

function Base.show(io::IO, m::GmshMesh)
    println(io, m.meshFormat)
    println(io, m.physicalnames)
    println(io, m.entities)
    println(io, m.nodeBlocks)
    println(io, m.elementblocks)
end






