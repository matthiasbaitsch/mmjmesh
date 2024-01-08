module Groups

using MMJMesh
import MMJMesh.MMJBase: SeqIntSet

struct EntityGroup <: AbstractVector{Int}
    name::Symbol
    pdim::Int
    indexes::SeqIntSet
end

end
