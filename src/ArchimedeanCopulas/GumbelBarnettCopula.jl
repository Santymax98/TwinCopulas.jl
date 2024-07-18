"""
    GumbelBarnettCopula{P}

Fields:

  - θ::Real - parameter

Constructor

    GumbelBarnettCopula(θ)

The bivariate Gumbel-Barnett copula is parameterized by ``\\theta \\in (0,1]``. It is an Archimedean copula with generator :

```math
\\phi(t) = \\exp{θ^{-1}(1-e^{t})}, 0 \\leq \\theta \\leq 1.
```

It has a few special cases: 
- When θ = 0, it is the IndependentCopula

References:
* Joe, H. (2014). Dependence modeling with copulas. CRC press, Page.437
* Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct GumbelBarnettCopula{T} <: ArchimedeanCopula{T}
    θ::T
    function GumbelBarnettCopula(θ)
        if !(0 <= θ <= 1)
            throw(ArgumentError("Theta must be in [0,1]"))
        elseif θ == 0
            return IndependentCopula()
        else
            return new{typeof(θ)}(θ)
        end 
    end
end

𝘙(C::GumbelBarnettCopula) = 1 # ϕ is not completely monotone so we cannot use the Williamsons transform

ϕ(C::GumbelBarnettCopula, x) = exp((1 - exp(x)) / C.θ)

ϕ⁻¹(C::GumbelBarnettCopula, x) = log1p(-C.θ * log(x))

dϕ(C::GumbelBarnettCopula, x) = -exp(x) / C.θ * exp((1 - exp(x)) / C.θ)

dϕ⁻¹(C::GumbelBarnettCopula, x) = -C.θ / (x * (1 - C.θ * log(x)))

d²ϕ(C::GumbelBarnettCopula, x) = -exp(x) / C.θ * exp((1 - exp(x)) / C.θ) + exp(2x) / C.θ^2 * exp((1 - exp(x)) / C.θ)
