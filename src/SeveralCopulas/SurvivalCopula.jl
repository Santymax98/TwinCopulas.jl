"""
    SurvivalCopula{P}

Fields:

    - C::bicopula
    
Constructor

    SurvivalCopula(C)

The bivariate Survival copula is defined as

```math
\\overbar{C}(u_1,u_2) = C(u_1, u_2) + u_1 + u_2 - 1
```
Where ``C(u_1, u_2)`` is a copula as we know them.

References:
* Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct SurvivalCopula{P<:bicopula} <: bicopula
    C::P
    function SurvivalCopula(C::P) where {P<:bicopula}
        new{P}(C)
    end
end

function Distributions.cdf(Cop::SurvivalCopula, x::AbstractVector) 
    u1, u2 = x
    copula = Distributions.cdf(Cop.C, [1-u1, 1-u2])
    return copula + u1 + u2 - 1
end

function Distributions.pdf(Cop::SurvivalCopula, x::AbstractVector)
    u1, u2 = x
    return Distributions.pdf(Cop.C, [1-u1, 1-u2])
end
function Distributions._rand!(rng::Distributions.AbstractRNG, Cop::SurvivalCopula, x::AbstractVector{T}) where {T<:Real}
    sample1, sample2 = Distributions.rand(rng, Cop.C, 2)
    x[1] = 1 - sample1
    x[2] = 1 - sample2
    return x
end