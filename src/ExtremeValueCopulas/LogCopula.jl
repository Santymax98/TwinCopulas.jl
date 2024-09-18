"""
    LogCopula{P}

Fields:

    - θ::Real - parameter
    
Constructor

    LogCopula(θ)

The bivariate Logistic copula (or Gumbel Copula) is parameterized by ``\\theta \\in [1,\\infty)``. It is an Extreme value copula with Pickands dependence function: 

```math
A(t) = (t^{\\theta}+(1-t)^{\\theta})^{\\frac{1}{\\theta}}
```

It has a few special cases: 
- When θ = 1, it is the IndependentCopula
- When θ = ∞, is is the MCopula (Upper Frechet-Hoeffding bound)

References:
* Bivariate extreme value theory: models and estimation. Biometrika, 1988.
"""
struct LogCopula{P} <: ExtremeValueCopula{P}
    θ::P  # Copula parameter
    function LogCopula(θ)
        if !(1 <= θ)
            throw(ArgumentError("El parámetro θ debe estar en [1, ∞)"))
        elseif θ == 1
            return IndependentCopula()
        elseif θ == Inf
            return MCopula()
        else
            return new{typeof(θ)}(θ)
        end
    end
end
# Definir la función ℓ específica para la copula de Galambos
function ℓ(G::LogCopula, t::Vector)
    θ = G.θ
    t₁, t₂ = t
    return (t₁^θ + t₂^θ)^(1/θ)
end
# Definir la función A específica para la copula de Galambos
function 𝘈(C::LogCopula, t::Real)
    θ = C.θ
    return (t^θ + (1 - t)^θ)^(1/θ)
end