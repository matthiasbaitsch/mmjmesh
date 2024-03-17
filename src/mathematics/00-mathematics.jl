module Mathematics


# Required modules
using IntervalSets
using IntervalSets
import Polynomials
using LinearAlgebra


# Parts
include("mappings.jl")
include("fefunctions.jl")


# Exports
export AbstractMapping, MappingFromR, FunctionToR, FunctionRToR
export domaintype, codomaintype, domain, valueat, derivativeat, derivative
export R, RPlus, R0Plus, IHat
export antiderivative, integrate, sample, plot, pois, roots
export Sin, Cos, Polynomial, fromroots, lagrangepolynomials, monomials, degree
export lagrangebasis1d, Interpolation
export ParametricCurve, UnitNormal


end