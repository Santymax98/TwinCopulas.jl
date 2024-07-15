"""
    WCopula

Constructor

    WCopula()

The [Lower Frechet-Hoeffding bound](https://en.wikipedia.org/wiki/Copula_(probability_theory)#Fr%C3%A9chet%E2%80%93Hoeffding_copula_bounds) is defined as

```math
W(u_1,u_2) = \\max\\{u_1 + u_2 - 1, 0\\} 
```
References:
* [nelsen2006](@cite) Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""

struct WCopula <: bicopula end

# Función CDF para la W Copula
function Distributions.cdf(C::WCopula, u::AbstractVector)
    return max(sum(u) - 1, 0)
end

# Función PDF para la W Copula
function Distributions.pdf(C::WCopula, u::AbstractVector)
    throw(ArgumentError("There is no density for the Frechet-Hoeffding bounds."))
end

function Distributions._rand!(rng::Distributions.AbstractRNG, C::WCopula, x::AbstractVector{T}) where {T<:Real}
    U = rand(rng, Distributions.Uniform(0,1))
    x[1] = U
    x[2] = 1-U
    return x
end