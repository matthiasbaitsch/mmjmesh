module Mathematics


# Required modules
using IntervalSets
using LinearAlgebra
using IntervalSets

import Polynomials
import CairoMakie as cm # XXX get rid of


# Parts
include("mappings.jl")


# Exports
export All
export AbstractMapping, MappingFromR, FunctionToR, FunctionRToR
export domaintype, codomaintype, domain, valueat, derivativeat, derivative
export R, RPlus, R0Plus, IHat
export antiderivative, integrate, sample, plot, pois, roots
export Sin, Cos, Polynomial, fromroots, lagrangepolynomials, monomials, degree


end