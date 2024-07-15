"""
    HuslerReissCopula{P}

Fields:

    - θ::Real - parameter
    
Constructor

    HuslerReissCopula(θ)

The bivariate Husler-Reiss copula is parameterized by ``\\theta \\in [0,\\infty)``. It is an Extreme value copula with Pickands dependence function: 

```math
A(t) = t\\Phi(\\theta^{-1}+\\frac{1}{2}\\theta\\log(\\frac{t}{1-t})) +(1-t)\\Phi(\\theta^{-1}+\\frac{1}{2}\\theta\\log(\\frac{1-t}{t}))
```
Where ``\\Phi``is the cumulative distribution function (CDF) of the standard normal distribution.

It has a few special cases:

- When θ = 0, it is the Independent Copula
- When θ = ∞, it is the M Copula (Upper Frechet-Hoeffding bound)

References:
* [Husler1989](@cite) Maxima of normal random vectors: between independence and complete dependence. Statist. Probab. 1989.
"""
struct HuslerReissCopula{P} <: ExtremeValueCopula{P}
    θ::P # Copula parameter
    function HuslerReissCopula(θ)
        if θ < 0
            throw(ArgumentError("Theta must be ≥ 0"))
        elseif θ == 0
            return IndependentCopula()
        elseif θ == Inf
            return MCopula()
        else
            return new{typeof(θ)}(θ)
        end
    end
end
# Definir la función ℓ específica para la copula de Galambos
function ℓ(H::HuslerReissCopula, t::Vector)
    θ = H.θ
    t₁, t₂ = t
    return t₁*Distributions.cdf(Distributions.Normal(),θ^(-1)+0.5*θ*log(t₁/t₂))+t₂*Distributions.cdf(Distributions.Normal(),θ^(-1)+0.5*θ*log(t₂/t₁))
end

# Definir la función A específica para la copula de Galambos
function 𝘈(H::HuslerReissCopula, t::Real)
    θ = H.θ
    term1 = t * Distributions.cdf(Normal(), θ^(-1) + 0.5 * θ * log(t / (1 - t)))
    term2 = (1 - t) * Distributions.cdf(Normal(), θ^(-1) + 0.5 * θ * log((1 - t) / t))
    
    A = term1 + term2
    
    
    return A
end

function d𝘈(H::HuslerReissCopula, t::Real)
    θ = H.θ
    # Derivada de A(x) respecto a t
    dA_term1 = Distributions.cdf(Normal(), θ^(-1) + 0.5 * θ * log(t / (1 - t))) + 
                  t * Distributions.pdf(Normal(), θ^(-1) + 0.5 * θ * log(t / (1 - t))) * (0.5 * θ * (1 / t + 1 / (1 - t)))
                  
    dA_term2 = -Distributions.cdf(Normal(), θ^(-1) + 0.5 * θ * log((1 - t) / t)) + 
                  (1 - t) * Distributions.pdf(Normal(), θ^(-1) + 0.5 * θ * log((1 - t) / t)) * (0.5 * θ * (-1 / t - 1 / (1 - t)))
    
    dA = dA_term1 + dA_term2
    
    return dA
end