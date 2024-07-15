"""
    JoeCopula{P}

Fields:

  - θ::Real - parameter

Constructor

    JoeCopula(θ)

The bivariate [Joe](https://en.wikipedia.org/wiki/Copula_(probability_theory)#Most_important_Archimedean_copulas) copula is parameterized by ``\\theta \\in [1,\\infty)``. It is an Archimedean copula with generator : 

```math
\\phi(t) = 1 - \\left(1 - e^{-t}\\right)^{\\frac{1}{\\theta}}
```

It has a few special cases: 
- When θ = 1, it is the IndependentCopula
- When θ = ∞, is is the MCopula (Upper Frechet-Hoeffding bound)

References:
* [nelsen2006](@cite) Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct JoeCopula{P} <: ArchimedeanCopula{P}
    θ::P
    function JoeCopula(θ)
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

𝘙(C::JoeCopula) = Sibuya(1 / C.θ)

ϕ(C::JoeCopula, x) = 1 - (-expm1(-x))^(1 / C.θ)

dϕ(C::JoeCopula, x) = -(1 / C.θ) * (-expm1(-x))^(1 / C.θ - 1) * exp(-x)

d²ϕ(C::JoeCopula, x) = -(1 / C.θ) * ((1 / C.θ - 1) * exp(-2 * x) * (-expm1(-x))^(1 / C.θ - 2) - exp(-x) * (-expm1(-x))^(1 / C.θ - 1))

ϕ⁻¹(C::JoeCopula, x) = -log1p(-(1 - x)^C.θ)

dϕ⁻¹(C::JoeCopula, x) = (C.θ * (1 - x)^(C.θ - 1)) / (1 - (1 - x)^C.θ)

τ(C::JoeCopula) = 1 - 4 * sum(1 / (k * (2 + k * C.θ) * (C.θ * (k - 1) + 2)) for k in 1:1000)

λᵤ(C::JoeCopula) = 2 - 2^(1 / C.θ)

λₗ(C::JoeCopula) = 0.0