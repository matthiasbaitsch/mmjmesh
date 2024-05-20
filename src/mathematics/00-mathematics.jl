module Mathematics


# Required modules
using DomainSets
using IntervalSets
using LinearAlgebra
import Polynomials


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
export AdHocMapping
export makefunction
export lagrangebasis1d, Interpolation
export ParametricCurve, UnitNormal
export ValueAtLF, DerivativeAtLF


end