module Mathematics

using Symbolics
using SymbolicUtils
using IntervalSets
using StaticArrays
using LinearAlgebra
using DomainSets: ×, Rectangle, ProductDomain, dimension, component
using SymbolicUtils: Postwalk, Chain
include("poly.jl") # My own version of this code - forgot why I did this

import Polynomials as P
import MMJMesh.Mathematics.FixedPolynomials as FP

using MMJMesh.MMJBase

include("mappings.jl")
include("forms.jl")
include("fefunctions.jl")
include("mpolynomials.jl")


# General concept of mapping
export MappingFromComponents
export AbstractMapping, MappingFromR, MappingFromRn, FunctionToR, FunctionRToR, FunctionRnToR
export domaintype, codomaintype, domain, valueat
export derivativeat, derivative, pderivativeat, pderivative
export R, R2, RPlus, R0Plus, IHat, ReferenceInterval, QHat, ReferenceQuadrilateral
export antiderivative, integrate, sample, plot, pois, roots
export Zero, One

# Functions Rn to R
export gradient, gradientat, hessian, hessianat, laplacian, laplacianat, ∇, H, Δ
export ProductFunction
export ∂x, ∂y, ∂xx, ∂yy, ∂xy

# Vector fields
export VectorField, div, divergence, divergenceat

# Parametric curves
export ParametricCurve, UnitNormal

# Special functions
export Sin, Cos
export MPolynomial, mmonomials, simplifyx
export AdHocMapping, makefunction
export AffineMapping
export Polynomial, Interpolation, fromroots, lagrangepolynomials, monomials, degree, lagrangebasis1d

# Forms
export ValueAtLF, DerivativeAtLF, DDerivativeAtLF, PDerivativeAtLF

end