"""
    ClaytonCopula{P}

Fields:

  - θ::Real - parameter

Constructor

    ClaytonCopula(d,θ)

The bivariate [Clayton](https://en.wikipedia.org/wiki/Copula_(probability_theory)#Most_important_Archimedean_copulas) copula is parameterized by ``\\theta \\in [-1,\\infty)``. It is an Archimedean copula with generator : 

```math
\\phi(t) = \\left(1+\\mathrm{sign}(\\theta)*t\\right)^{-1\\frac{1}{\\theta}}
```

It has a few special cases: 
- When θ = -1, it is the WCopula (Lower Frechet-Hoeffding bound)
- When θ = 0, it is the IndependentCopula
- When θ = 1, it is the UtilCopula
- When θ = ∞, is is the MCopula (Upper Frechet-Hoeffding bound)

References:
* [nelsen2006](@cite) Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct ClaytonCopula{P} <: ArchimedeanCopula{P}
    θ::P
    function ClaytonCopula(θ) 
        if !(-1 <= θ) 
            throw(ArgumentError("El parámetro θ debe estar en [-1, ∞)")) 
        elseif θ == -1 
            return WCopula() 
        elseif θ == 0 
            return IndependentCopula() 
        elseif θ == 1 
            return UtilCopula() 
        elseif θ == Inf 
            return MCopula() 
        else 
            return new{typeof(θ)}(θ) 
        end 
    end
end

ϕ(C::ClaytonCopula, x) = max(1 + C.θ * x, zero(x))^(-1 / C.θ)

dϕ(C::ClaytonCopula, x) = (-1 / C.θ) * (1 + x)^(-1 / C.θ - 1)

d²ϕ(C::ClaytonCopula, x) = (1 / C.θ) * (1 / C.θ + 1) * (1 + x)^(-1 / C.θ - 2)

ϕ⁻¹(C::ClaytonCopula, x) = (x^(-C.θ) - 1) / C.θ

dϕ⁻¹(C::ClaytonCopula, x) = -x^(-C.θ - 1)

τ(C::ClaytonCopula) = C.θ / (C.θ + 2)

λᵤ(C::ClaytonCopula) = 0.0

λₗ(C::ClaytonCopula) = 2^(-1 / C.θ)

𝘙(C::ClaytonCopula) = C.θ > 0 ? Distributions.Gamma(1/C.θ) : 1
