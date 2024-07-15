"""
    MardiaCopula{P}

Fields:

    - θ::Real - parameter
    
Constructor

    MardiaCopula(θ)

The bivariate Mardia copula is parameterized by ``\\theta \\in [-1,1]``. It is constructed as: 

```math
C(u_1, u_2) = \\frac{\\theta^2(1+\\theta)}{2}\\min\\{u_1,u_2\\} + (1-\\theta^2)u_1u_2 + \\frac{\\theta^2(1-\\theta)}{2}\\max\\{u_1+u_2-1,0\\}
```

It has a few special cases: 
- When θ = 0, it is the Independent Copula
- When θ = 1, it is the MCopula (Upper Frechet-Hoeffding bound) 
- When θ = -1, is is the WCopula (Lower Frechet-Hoeffding bound)

References:
* [nelsen2006](@cite) Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct MardiaCopula{P} <: bicopula
    θ::P  # Copula parameter
    function MardiaCopula(θ)
        if !(-1 <= θ <= 1)
            throw(ArgumentError("Theta must be in [-1,1]"))
        elseif θ == 0
            return IndependentCopula()
        elseif θ == 1
            return MCopula()
        elseif θ == -1
            return WCopula()
        else
            return new{typeof(θ)}(θ)
        end
    end
end

function Distributions.cdf(C::MardiaCopula, u::AbstractVector)
    θ = C.θ
    u1, u2 = u
    term1 = (θ^2 * (1 + θ) / 2) * min(u1, u2)
    term2 = (1 - θ^2) * u1 * u2
    term3 = (θ^2 * (1 - θ) / 2) * max(u1 + u2 - 1, 0)
    return term1 + term2 + term3
end

function Distributions.pdf(C::MardiaCopula, u::AbstractVector)
    throw(ArgumentError("There is no density for the Frechet Copula."))
end

function Distributions._rand!(rng::Distributions.AbstractRNG, C::MardiaCopula, x::AbstractVector{T}) where {T<:Real}
    θ = C.θ
    u1, u2 = rand(rng, Distributions.Uniform(0,1), 2)
    p = [θ^2 * (1 + θ) / 2, 1 - θ^2, θ^2 * (1 - θ) / 2]
    z = rand(rng, Distributions.Categorical(p))
    
    if z == 1
        u1 = min(u1, u2)
        x[1] = u1
        x[2] = u1
    elseif z == 2
        x[1] = u1
        x[2] = u2
    else
        u1 = max(u1 + u2 - 1, 0)
        x[1] = u1
        x[2] = 1 - u1
    end
    
    return x
end