"""
    GumbelCopula{P}

Fields:

  - θ::Real - parameter

Constructor

    GumbelCopula(d,θ)

The bivariate [Gumbel](https://en.wikipedia.org/wiki/Copula_(probability_theory)#Most_important_Archimedean_copulas) copula is parameterized by ``\\theta \\in [1,\\infty)``. It is an Archimedean copula with generator : 

```math
\\phi(t) = \\exp{-t^{\\frac{1}{θ}}}
```

It has a few special cases: 
- When θ = 1, it is the IndependentCopula
- When θ = ∞, is is the MCopula (Upper Frechet-Hoeffding bound)

References:
* Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct GumbelCopula{P} <: ArchimedeanCopula{P}
    θ::P
    function GumbelCopula(θ) 
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

𝘙(C::GumbelCopula) = Stable(1 / C.θ)

ϕ(C::GumbelCopula, x) = exp(-x^(1 / C.θ))

dϕ(C::GumbelCopula, x) = - (1 / C.θ) * x^(1 / C.θ - 1) * exp(-x^(1 / C.θ))

d²ϕ(C::GumbelCopula, x) = -(1 / C.θ) * (1 / C.θ - 1) * x^(1 / C.θ - 2) * exp(-x^(1 / C.θ)) + (1 / C.θ^2) * x^(2 / C.θ - 2) * exp(-x^(1 / C.θ))

ϕ⁻¹(C::GumbelCopula, x) = (-log(x))^C.θ

dϕ⁻¹(C::GumbelCopula, x) = -(C.θ / x) * (-log(x))^(C.θ - 1)

τ(C::GumbelCopula) = (C.θ - 1) / C.θ

λᵤ(C::GumbelCopula) = 2 - 2^(1 / C.θ)

λₗ(C::GumbelCopula) = 0.0
