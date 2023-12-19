module Mathematics

# Modules needed by this module
using IntervalSets

# Functions extended by this module
import Base.in
import Base.isequal
import Base.isapprox

using IntervalSets
import Polynomials
import CairoMakie as cm


# Parts
include("mappings.jl")
export Everything, AbstractMapping, AbstractFunctionRToR, Sin, Cos, Polynomial, monomials, plot

# Exports
## regions.jl
# export Region, Interval, dim, discretize

end