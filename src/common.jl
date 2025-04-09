# -------------------------------------------------------------------------------------------------
# Operations for various mesh components
# -------------------------------------------------------------------------------------------------


"""
    id(o)

Returns the numerical id of an object
"""
function id end

"""
    pdim(o)

Parametric dimension of object `o`.
"""
function pdim end

"""
    gdim(o)

Dimension of space object `o` is embedded in.
"""
function gdim end

"""
    coordinate(o, i)

Returns the `i`th coordinate of `o`.
"""
function coordinate end

"""
    coordinates(o)
    coordinates(o, s)
    coordinates(o, s, c)

Returns coordinates of `o`.
"""
function coordinates end

"""
    dimension(o)

The dimension of object `o`. 
"""
function dimension end

"""
    mesh(o)

Returns the mesh, object `o` belongs to.
"""
function mesh end


## Export

export id
export pdim, gdim
export coordinate, coordinates
export dimension
export mesh


# -------------------------------------------------------------------------------------------------
# Common vector types
# -------------------------------------------------------------------------------------------------

const RealVec = AbstractVector{<:Real}
const RealMat = AbstractMatrix{<:Real}
const RealVecOrMat = Union{RealVec,RealMat}
const RealVecVec = AbstractVector{<:AbstractVector{<:Real}}
const IntegerVec = AbstractVector{<:Integer}
const IntegerMat = AbstractMatrix{<:Integer}
const IntegerVecOrMat = Union{IntegerVec,IntegerMat}
const IntegerVecVec = AbstractVector{<:AbstractVector{<:Integer}}


## Export

export RealVec, RealMat, RealVecOrMat, RealVecVec
export IntegerVec, IntegerMat, IntegerVecOrMat, IntegerVecVec
