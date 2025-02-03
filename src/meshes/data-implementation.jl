
function setdata!(m::Mesh, id::Symbol, value)
    m.data.meshdata[id] = value
    nothing
end

function setdata!(g::Group, id::Symbol, value)
    g.mesh.data.groupdata[Pair(id, name(g))] = value
    nothing
end

function setdata!(e::MeshEntity, id::Symbol, value)
    e.mesh.data.entitydata[(id, pdim(e), index(e))] = value
    nothing
end

data(m::Mesh, id::Symbol) = m.data.meshdata[id]

data(g::Group, id::Symbol) = g.mesh.data.groupdata[Pair(id, name(g))]

hasdata(o, id::Symbol) = !isnothing(data(o, id))


_retrieve_entitydata(::MeshEntity, v::Number) = v

function _retrieve_entitydata(e::MeshEntity, v::AbstractVector)
    @assert nentities(e.mesh, pdim(e)) == length(v)
    return v[index(e)]
end


function data(e::MeshEntity, id::Symbol)
    d = e.mesh.data

    # Entity
    key = (id, pdim(e), index(e))
    haskey(d.entitydata, key) && return _retrieve_entitydata(e, d.entitydata[key])

    # Group
    for gn in groupnames(e)
        key = Pair(id, gn)
        haskey(d.groupdata, key) && return _retrieve_entitydata(e, d.groupdata[key])
    end

    # Mesh
    haskey(d.meshdata, id) && return _retrieve_entitydata(e, d.meshdata[id])

    nothing
end