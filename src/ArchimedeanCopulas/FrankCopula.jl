"""
    FrankCopula{P}

Fields:

  - θ::Real - parameter

Constructor

    FrankCopula(θ)

The bivariate [Frank](https://en.wikipedia.org/wiki/Copula_(probability_theory)#Most_important_Archimedean_copulas) copula is parameterized by ``\\theta \\in (-\\infty,\\infty)``. It is an Archimedean copula with generator : 

```math
\\phi(t) = -\\frac{\\log\\left(1+e^{-t}(e^{-\\theta-1})\\right)}{\theta}
```

It has a few special cases: 
- When θ = -∞, it is the WCopula (Lower Frechet-Hoeffding bound)
- When θ = 1, it is the IndependentCopula
- When θ = ∞, is is the MCopula (Upper Frechet-Hoeffding bound)

References:
* [nelsen2006](@cite) Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct FrankCopula{P} <: ArchimedeanCopula{P}
    θ::P
    function FrankCopula(θ) 
        if θ == -Inf 
            return WCopula() 
        elseif θ == 0 
            return IndependentCopula() 
        elseif θ == Inf 
            return MCopula() 
        else 
            return new{typeof(θ)}(θ) 
        end 
    end
end
𝘙(C::FrankCopula) = C.θ > 0 ? Logarithmic(1 - exp(-C.θ)) : 1

ϕ(C::FrankCopula, x) = C.θ > 0 ? -LogExpFunctions.log1mexp(LogExpFunctions.log1mexp(-C.θ) - x) / C.θ : -log1p(exp(-x) * expm1(-C.θ)) / C.θ

dϕ(C::FrankCopula, x) = (exp(-x) * expm1(-C.θ)) / (C.θ * (exp(-x) * expm1(-C.θ)) + 1)

d²ϕ(C::FrankCopula, x) = -(exp(-x) * expm1(-C.θ)) / (C.θ * (exp(-x) * expm1(-C.θ) + 1)^2)

ϕ⁻¹(C::FrankCopula, x) = -log((exp(-C.θ * x) - 1) / (exp(-C.θ) - 1))

dϕ⁻¹(C::FrankCopula, x) = (C.θ * exp(-C.θ * x)) / (exp(-C.θ * x) - 1)

τ(C::FrankCopula) = (1 + 4 * (𝘋(1, C.θ) - 1) / C.θ)

ρₛ(C::FrankCopula) = (1 + 12 * (𝘋(2, C.θ) - 𝘋(1, C.θ)) / C.θ)

λᵤ(C::FrankCopula) = 0.0

λₗ(C::FrankCopula) = 0.0

function 𝘋(n, x)
    if x == 0
        return 1.0
    else
        integrand(t) = (t^n) / (exp(t) - 1)
        integral, _ = QuadGK.quadgk(integrand, 0, x)
        return (n / x^n) * integral
    end
end