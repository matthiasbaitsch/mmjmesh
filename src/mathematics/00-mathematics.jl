module Mathematics

using Symbolics
using SymbolicUtils
using IntervalSets
using StaticArrays
using LinearAlgebra
using SymbolicUtils: Postwalk, Chain
using DomainSets: ×, ProductDomain, Rectangle

import DomainSets
import Polynomials

# Multivariate polynomials, my own version of 
## FixedPolynomials - forgot why I did this
include("poly.jl")
import MMJMesh.Mathematics.FixedPolynomials as FP

using MMJMesh.MMJBase

import MMJMesh: dimension

include("domains.jl")
include("mappings.jl")
include("forms.jl")
include("mpolynomials.jl")
include("interpolations.jl")
include("parametertypes.jl")

include("spaces.jl")
include("finiteelements.jl")

# Domains
export R, RPlus, R⁺, R0Plus, R⁺₀, IHat, ReferenceInterval
export R2, R², R3, R³, QHat, ReferenceQuadrilateral
export InR, InRⁿ, InR2, InR², InR3, InR³, InRnxm, InRⁿˣᵐ
export points, dimension

# General concept of mapping
export AbstractMapping, MappingFromR, MappingFromRn, FunctionToR
export FunctionRToR, FunctionRnToR, MappingRnToRm
export domaintype, codomaintype, domain, degree, degrees, valueat
export derivativeat, derivative, pderivativeat, pderivative
export antiderivative, integrate, sample, plot, pois, roots

# Composed functions
export MappingFromComponents

# Special functions I
export Zero, One

# Functions Rn to R
export gradient, gradientat, hessian, hessianat, laplacian, laplacianat, ∇, H, Δ
export ProductFunction
export ∂x, ∂y, ∂xx, ∂yy, ∂xy

# Mappings Rn to Rm
export jacobian, jacobianat

# Vector fields
export VectorField, div, divergence, divergenceat

# Parametric curves
export ParametricCurve, UnitNormal

# Special functions II
export Sin, Cos
export MPolynomial, mmonomials, simplifyx
export AdHocMapping, makefunction
export AffineMapping
export Polynomial, Interpolation, fromroots, lagrangepolynomials, monomials, lagrangebasis1d

# Forms
export Form, LinearForm, BilinearForm, ValueAtLF, DerivativeAtLF, DDerivativeAtLF, PDerivativeAtLF
export ∂xLF, ∂yLF, ∂xyLF

# Spaces
export FunctionSpace, PolynomialSpace, P, Q, S, Q23R, basis

# Finite elements
export FiniteElement, makeelement, nodalbasis

end