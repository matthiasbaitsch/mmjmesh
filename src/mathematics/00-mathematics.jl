module Mathematics

# Modules needed by this module
using IntervalSets
using LinearAlgebra

# Functions extended by this module
import Base.in
import Base.isequal
import Base.isapprox

using IntervalSets
import Polynomials
import CairoMakie as cm


# Parts
include("mappings.jl")

# XXX
export Everything, AbstractMapping, AbstractFunctionRToR, Sin, Cos, Polynomial, monomials, plot

export derivative, antiderivative, integrate, fromroots, lagrangepolynomials

# Exports
## regions.jl
# export Region, Interval, dim, discretize

end