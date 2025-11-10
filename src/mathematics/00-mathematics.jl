module Mathematics

using QuadGK
using Symbolics
using IntervalSets
using StaticArrays
using SymbolicUtils
using LinearAlgebra
using SymbolicUtils: Postwalk, Chain
using DomainSets: ×, ProductDomain, Rectangle

import DomainSets
import Polynomials

using MMJMesh
using MMJMesh.MMJBase

include("domains.jl")
include("mappings.jl")
include("forms.jl")
include("mpolynomials.jl")
include("curves.jl")
include("interpolations.jl")
include("parametertypes.jl")
include("piecewise.jl")
include("spaces.jl")
include("finiteelements.jl")

# Domains
export R, RPlus, R⁺, R0Plus, R⁺₀, IHat, ReferenceInterval
export R2, R², R3, R³, QHat, ReferenceQuadrilateral
export InR, InRⁿ, InR2, InR², InR3, InR³, InRmxn, InRᵐˣⁿ
export points, dimension

# General concept of mapping
export AbstractMapping, MappingFromR, MappingFromRn
export FunctionToR, MappingToRn, MappingToRmxn
export FunctionRToR, FunctionRnToR, MappingRnToRm
export domaintype, codomaintype, domain, degree, degrees
export constval, valueat
export derivativeat, derivative, pderivativeat, pderivative
export antiderivative, integrate, sample, plot, pois, roots

# Composed functions
export MappingFromComponents

# Special functions I
export Zero, One, Identity

# Variable
export parameter

# Functions Rn to R
export gradient, gradientat, hessian, hessianat, laplacian, laplacianat, ∇, H, Δ
export ProductFunction
export ∂x, ∂y, ∂xx, ∂yy, ∂xy

# Mappings to Rn
export components

# Mappings Rn to Rm
export jacobian, jacobianat

# Vector fields
export VectorField, div, divergence, divergenceat

# Parametric curves
export ParametricCurve, UnitNormal

# Special functions II
export constant
export Sin, Cos, Exp
export mmonomials, affinefunction, simplifyx
export AdHocMapping, makefunction
export AffineMapping
export Polynomial, Interpolation
export monomial, monomials, fromroots, lagrangepolynomials, lagrangebasis1d
export PiecewiseFunction, npieces, interpolate

# Polynomials of multiple variables
export MPolynomial, PolynomialRnToR, PolynomialRnToRm, PolynomialRnToRpxq, 
    nterms, exponents, coefficient, coefficients, degrees, mmonomials

# Forms
export Form, LinearForm, BilinearForm, ValueAtLF, DerivativeAtLF, DDerivativeAtLF, PDerivativeAtLF
export ∂xLF, ∂yLF, ∂xyLF

# Spaces
export FunctionSpace, PolynomialSpace, P, Q, S, Q23R, basis

# Curves
export linesegment, polynomialinterpolation

# Finite elements
export hatfunctions, FiniteElement, makeelement, nodalbasis

end