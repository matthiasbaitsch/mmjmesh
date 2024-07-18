module Mathematics

using Symbolics
using SymbolicUtils
using IntervalSets
using StaticArrays
using LinearAlgebra
using SymbolicUtils: Postwalk, Chain
using DomainSets: ×, ProductDomain, Rectangle

# Multivariate polynomials, my own version of 
## FixedPolynomials - forgot why I did this
include("poly.jl")

import DomainSets
import Polynomials
import MMJMesh.Mathematics.FixedPolynomials as FP

using MMJMesh.MMJBase

include("domains.jl")
include("mappings.jl")
include("forms.jl")
include("fefunctions.jl")
include("mpolynomials.jl")
include("parametertypes.jl")

# Domains
export R, RPlus, R⁺, R0Plus, R⁺₀, IHat, ReferenceInterval
export R2, R², R3, R³, QHat, ReferenceQuadrilateral
export InR, InRⁿ, InR2, InR², InR3, InR³, InRⁿˣᵐ
export dimension

# General concept of mapping
export MappingFromComponents
export AbstractMapping, MappingFromR, MappingFromRn, FunctionToR, FunctionRToR, FunctionRnToR
export domaintype, codomaintype, domain, valueat
export derivativeat, derivative, pderivativeat, pderivative
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
export Polynomial, Interpolation, fromroots, lagrangepolynomials, monomials, lagrangebasis1d

export degree, degrees

# Forms
export ValueAtLF, DerivativeAtLF, DDerivativeAtLF, PDerivativeAtLF

end