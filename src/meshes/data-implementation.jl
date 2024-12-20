
function setdata!(m::Mesh, id::Symbol, value)
    m.data.meshdata[id] = value
end

function setdata!(g::Group, id::Symbol, value)
    g.mesh.data.groupdata[Pair(id, name(g))] = value
end

function setdata!(e::MeshEntity, id::Symbol, value)
    e.mesh.data.entitydata[(id, pdim(e), index(e))] = value
end

data(m::Mesh, id::Symbol) = m.data.meshdata[id]

data(g::Group, id::Symbol) = g.mesh.data.groupdata[Pair(id, name(g))]

hasdata(o, id::Symbol) = !isnothing(data(o, id))

function data(e::MeshEntity, id::Symbol)
    d = e.mesh.data

    # Entity
    key = (id, pdim(e), index(e))
    haskey(d.entitydata, key) && return d.entitydata[key]

    # Group
    for gn in groupnames(e)
        key = Pair(id, gn)
        haskey(d.groupdata, key) && return d.groupdata[key]
    end

    # Mesh
    haskey(d.meshdata, id) && return d.meshdata[id]

    return nothing
end