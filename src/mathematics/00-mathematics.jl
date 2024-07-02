module Mathematics


include("poly.jl")
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
export Zero

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
export Polynomial, Interpolation, fromroots, lagrangepolynomials, monomials, degree, lagrangebasis1d

# Forms
export ValueAtLF, DerivativeAtLF, DDerivativeAtLF, PDerivativeAtLF

end