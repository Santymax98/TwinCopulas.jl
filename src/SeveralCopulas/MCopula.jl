"""
    MCopula

Constructor

    MCopula()

The [Upper Frechet-Hoeffding bound](https://en.wikipedia.org/wiki/Copula_(probability_theory)#Fr%C3%A9chet%E2%80%93Hoeffding_copula_bounds) is defined as

```math
M(u_1,u_2) = \\min\\{u_1,u_2\\}
```

References:
* [nelsen2006](@cite) Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct MCopula <: bicopula end

# Función CDF para la W Copula
function Distributions.cdf(C::MCopula, u::AbstractVector)
    u1, u2 = u
    return min(u1,u2)
end

# Función PDF para la W Copula
function Distributions.pdf(C::MCopula, u::AbstractVector)
    throw(ArgumentError("There is no density for the Frechet-Hoeffding bounds."))
end

function Distributions._rand!(rng::Distributions.AbstractRNG, C::MCopula, x::AbstractVector{T}) where {T<:Real}
    U = rand(rng, Distributions.Uniform(0,1))
    x[1] = U
    x[2] = U
    return x
end