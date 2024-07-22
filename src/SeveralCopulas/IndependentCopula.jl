"""
    IndependentCopula

Constructor

    IndependentCopula()

The bivariate [Independent Copula](https://en.wikipedia.org/wiki/Copula_(probability_theory)#Most_important_Archimedean_copulas) is
the simplest copula, that has the form : 

```math
C(u_1,u_2) = u_1u_2
```

References:
* Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct IndependentCopula <: bicopula end

# Función CDF para la Copula Independiente
function Distributions.cdf(C::IndependentCopula, u::AbstractArray)
    return prod(u)
end

# Función PDF para la Copula Independiente
function Distributions.pdf(C::IndependentCopula, u::AbstractVector)
    return 1.0
end

function Distributions._rand!(rng::Distributions.AbstractRNG, C::IndependentCopula, x::AbstractVector{T}) where {T<:Real}
    u1, u2 = rand(rng, Distributions.Uniform(0,1),2)
    x[1] = u1
    x[2] = u2
    return x
end