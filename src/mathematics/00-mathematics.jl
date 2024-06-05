module Mathematics


# Required modules
using DomainSets
using IntervalSets
using LinearAlgebra


# Parts
include("poly.jl")
include("mappings.jl")
include("forms.jl")
include("fefunctions.jl")


# Exports
export MappingFromComponents
export AbstractMapping, MappingFromR, MappingFromRn, FunctionToR, FunctionRToR, FunctionRnToR
export domaintype, codomaintype, domain, valueat
export derivativeat, derivative, pderivativeat, pderivative
export R, R2, RPlus, R0Plus, IHat, ReferenceInterval, QHat, ReferenceQuadrilateral
export antiderivative, integrate, sample, plot, pois, roots
export Sin, Cos, Polynomial, fromroots, lagrangepolynomials, monomials, degree

export AdHocMapping
export makefunction
export lagrangebasis1d, Interpolation

# Functions Rn to R
export gradient, gradientat, hessian, hessianat, laplacian, laplacianat, ∇, H, Δ
export ProductFunction, MPolynomial, mmonomials

# Vector fields
export VectorField, div, divergence, divergenceat

# Parametric curves
export ParametricCurve, UnitNormal

# Forms
export ValueAtLF, DerivativeAtLF, DDerivativeAtLF, PDerivativeAtLF


end