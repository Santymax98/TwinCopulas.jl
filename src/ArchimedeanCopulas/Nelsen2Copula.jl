"""
    Nelsen2Copula{P}

Fields:
  - θ::Real - parameter

Constructor

    Nelsen2Copula(θ)

The bivariate Nelsen2Copula copula is parameterized by ``\\theta \\in [1,\\infty)``. It is an Archimedean copula with generator : 

```math
\\phi(t) = 1 - t^{\\frac{1}{\\theta}}
```

It has a few special cases: 
- When θ = 1, it is the IndependentCopula
- When θ = ∞, is is the MCopula (Upper Frechet-Hoeffding bound)

References:
* Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct Nelsen2Copula{P} <: ArchimedeanCopula{P}
    θ::P
    function Nelsen2Copula(θ)
        if !(1 <= θ)
            throw(ArgumentError("El parámetro θ debe estar en [1, ∞]"))
        elseif θ == 1
            return WCopula()
        elseif θ == Inf
            return MCopula()
        else
            return new{typeof(θ)}(θ)
        end
    end
end

𝘙(C::Nelsen2Copula) = RadialDist(C)
ϕ(C::Nelsen2Copula, x) = max(1 - x^(1/C.θ), zero(x))
ϕ⁻¹(C::Nelsen2Copula, x) = (1 - x)^C.θ
dϕ(C::Nelsen2Copula, x) = -(1 / C.θ) * x^(1/C.θ - 1)
dϕ⁻¹(C::Nelsen2Copula, x) = C.θ * (1 - x)^(C.θ - 1)
d²ϕ(C::Nelsen2Copula, x) = -(1 / C.θ)*(1/C.θ - 1)*(x^(1/C.θ - 2))

τ(C::Nelsen2Copula) = (3 * C.θ - 2) / (3 * C.θ) - (2 * (1 - C.θ)^2 / (3 * C.θ^2)) * log(1 - C.θ)

ρₛ(C::Nelsen2Copula) = (12 * (1 + C.θ) / C.θ^2) * PolyLog.reli2(1 - C.θ) - (24 * (1 - C.θ) / C.θ^2) * log(1 - C.θ) - (3 * (C.θ + 12) / C.θ)

λᵤ(C::Nelsen2Copula) = 0.0

λₗ(C::Nelsen2Copula) = 0.0