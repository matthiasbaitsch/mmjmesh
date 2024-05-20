module Mathematics


# Required modules
using IntervalSets
using IntervalSets
import Polynomials
using LinearAlgebra


# Parts
include("mappings.jl")
include("forms.jl")
include("fefunctions.jl")


# Exports
export AbstractMapping, MappingFromR, MappingFromRn, FunctionToR, FunctionRToR, FunctionRnToR
export domaintype, codomaintype, domain, valueat, derivativeat, derivative
export R, RPlus, R0Plus, IHat
export antiderivative, integrate, sample, plot, pois, roots
export Sin, Cos, Polynomial, fromroots, lagrangepolynomials, monomials, degree
export lagrangebasis1d, Interpolation
export ParametricCurve, UnitNormal
export ValueAtLF, DerivativeAtLF


end