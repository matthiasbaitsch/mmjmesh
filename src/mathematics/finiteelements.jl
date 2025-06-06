"""
    hatfunctions(x::AbstractVector{<:Real})

Given an ascending node vector ``x_1, x_2, \\dots, x_{N_N}``, returns the 1D piecewise 
affine finite element basis functions ``\\varphi_1, \\dots, \\varphi_{N_N}`` with
``\\varphi_i(x_j) = \\delta_{ij}`` for ``i,j = 1, \\dots, N_N``.
"""
function hatfunctions(x::AbstractVector{<:Real})
    N = length(x)
    ei(i) = [i == j for j = 1:N]
    return [interpolate(x, ei(i), order=1) for i = 1:N]
end


"""
    FiniteElement(K, P, N)

Construct a finite element ``(K, \\mathcal{P}, \\mathcal{N})``  according to
Ciarlet's definition. Parameters are

- the domain `K`,

- the space `P` of functions defined on `K`,

- a vector of linear forms `N`,

- a `cache` holding computed objects such as a nodal basis.

"""
struct FiniteElement
    K
    P::Type{<:FunctionSpace}
    N::AbstractVector{<:LinearForm}
    cache::Dict{Symbol,Any}

    function FiniteElement(K, P::Type{<:FunctionSpace}, N::Vector{<:LinearForm})
        @assert dimension(P) == length(N)
        return new(K, P, N, Dict{Symbol,Any}())
    end
end


"""
    nodalbasis(e)
    nodalbasis(t, d)

Generate a nodal basis for element `e` defined on `d`. This default implementation 
might be overridden by specific element types. Option to specify domain only useful 
for symbolic domains.
"""
function nodalbasis(e::FiniteElement)
    if :nodalbasis ∉ keys(e.cache)
        # Hack to handle domains with symbolic limits
        d = e.K

        if !(d |> eltype |> eltype |> isbitstype)
            d = R^dimension(e.K)
        end

        # Generate nodal basis
        ps = basis(e.P, d)
        M = [n(p) for p in ps, n in e.N]
        e.cache[:nodalbasis] = inv(M) * ps
    end
    return e.cache[:nodalbasis]
end

nodalbasis(type::Symbol, domain; k = -1, a...) = makeelement(type, domain; k, a...) |> nodalbasis


"""
    _combineforms(points, genforms)

Helper function to construct linear forms.
"""
function _combineforms(points, genforms)
    if points isa Vector{<:Vector}
        points = reduce(vcat, points)
    end
    if genforms isa DataType
        genforms = [genforms]
    end
    return vec([gf(p) for gf in genforms, p in points])
end


# -------------------------------------------------------------------------------------------------
# List of elements
# -------------------------------------------------------------------------------------------------

struct ElementDescriptor
    description::String
    domaintypes::Vector{Type}
    degree::AbstractInterval{<:Integer}
    namedparameters::Dict{Symbol,Any}
    makeelement::Function
end

const elementdescriptors = Dict{Symbol,ElementDescriptor}()
const elementcache = Dict{Tuple{Symbol,Any,Integer,Any},FiniteElement}()
const elementdomainstocache = Set([IHat, QHat])

function addelementtype!(
    typeid; description, domaintypes, makeelement,
    degree=(1 .. 10000), namedparameters=Dict{Symbol,Any}()
)
    elementdescriptors[typeid] =
        ElementDescriptor(
            description,
            domaintypes,
            degree,
            namedparameters,
            makeelement
        )
end

function hasdomaintype(type, types)
    for t in types
        if type isa t
            return true
        end
    end
    return false
end

function domakeelement(K, d::ElementDescriptor, k::Integer, a...)
    @assert hasdomaintype(K, d.domaintypes)
    if k == -1
        return d.makeelement(K; a...)
    else
        @assert k ∈ d.degree
        return d.makeelement(K, k; a...)
    end
end

function makeelement(id::Symbol, K; k=-1, a...)
    d = elementdescriptors[id]

    if K ∈ elementdomainstocache
        key = (id, K, k, a)
        if key ∉ keys(elementcache)
            elementcache[key] = domakeelement(K, d, k, a...)
        end
        return elementcache[key]
    else
        return domakeelement(K, d, k, a...)
    end
end


# -------------------------------------------------------------------------------------------------
# Concrete element formulations
# -------------------------------------------------------------------------------------------------

"""
    lagrangeelement(K, k)

Lagrange type finite element.
"""
function lagrangeelement(K, k::Integer)
    @assert k >= 1
    return FiniteElement(
        K,
        Q{dimension(K),k},
        _combineforms(
            [
                points(K, :corners),
                points(K, :sides, k - 1),
                points(K, :interior, k - 1)
            ],
            ValueAtLF
        )
    )
end

addelementtype!(
    :lagrange, description="Standard Lagrange element",
    domaintypes=[Interval, Rectangle],
    makeelement=lagrangeelement
)

# Put in Cache
makeelement(:lagrange, IHat, k=1);
makeelement(:lagrange, QHat, k=1);


"""
    serendipityelement(K, k)

Serendipity element (currently only for rectangles and ``k=1,2``).
"""
function serendipityelement(K::Rectangle, k::Integer)
    @assert 1 <= k <= 2
    @assert dimension(K) == 2
    return FiniteElement(
        K,
        S{2,k},
        _combineforms(
            [
                points(K, :corners),
                points(K, :sides, k - 1)
            ],
            ValueAtLF
        )
    )
end

addelementtype!(
    :serendipity, description="Serendipity element",
    domaintypes=[Rectangle],
    makeelement=serendipityelement,
    degree=1 .. 2
)



"""
    hermiteelement(K; conforming=true)

Hermite type finite element. Generates the Bogner-Fox-Schmit element (Braess, p. 72) if `K` is a rectangle 
and `conforming` is true.
"""
hermiteelement(K::Interval) = FiniteElement(
    K,
    Q{1,3},
    _combineforms(points(K, :corners), [ValueAtLF, DerivativeAtLF])
)

function hermiteelement(K::Rectangle; conforming=true)
    if conforming
        return FiniteElement(
            K,
            Q{2,3},
            _combineforms(points(K, :corners), [ValueAtLF, ∂xLF, ∂yLF, ∂xyLF])
        )
    else
        return FiniteElement(
            K,
            Q23R,
            _combineforms(points(K, :corners), [ValueAtLF, ∂xLF, ∂yLF])
        )
    end
end

addelementtype!(
    :hermite,
    description="Hermite element",
    domaintypes=[Interval, Rectangle],
    makeelement=hermiteelement
)

# Put in Cache
makeelement(:hermite, IHat);
makeelement(:hermite, QHat);

