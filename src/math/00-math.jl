module Math

# Modules needed by this module

# Functions extended by this module
import Base.in
import Base.isequal
import Base.isapprox

# Parts
include("regions.jl")

# Exports
## regions.jl
export Region, Interval, dim, discretize

end