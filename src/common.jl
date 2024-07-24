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
    coordinates(os, i)

Returns the coordinates of `o` or the coordinates of the `i`th 
object in the collection `os`.
"""
function coordinates end

"""
    dimension(o)

The dimension of object `o`. 
"""
function dimension end


export pdim, gdim
export coordinate, coordinates
export dimension
