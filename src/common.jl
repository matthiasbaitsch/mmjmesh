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

TODO
"""
function coordinate end

"""
    coordinates(o)
    coordinates(o, i)

TODO
"""
function coordinates end

export pdim, gdim
export coordinate, coordinates
