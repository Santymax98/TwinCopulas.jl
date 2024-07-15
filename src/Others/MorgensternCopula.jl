"""
    MorgensternCopula{P}

Fields:

    - θ::Real - parameter
    
Constructor

    MorgensternCopula(θ)

The bivariate Morgenstern copula or FGM Copula is parameterized by ``\\theta \\in [-1,1]``. It is constructed as: 

```math
C(u_1, u_2) = u_1u_2(1+\\theta(1-u_1)(1-u_2))
```
It has a few special cases: 
- When θ = 0, it is the Independent Copula

References:
* [Joe1997](@cite) Joe, Harry, Multivariate Models and Multivariate Dependence Concepts, Chapman & Hall. 1997.
"""
struct MorgensternCopula{P} <: bicopula
    θ::P  # Copula parameter
    function MorgensternCopula(θ)
        if !(-1 <= θ <= 1)
            throw(ArgumentError("Theta must be in [-1,1]"))
        elseif θ == 0
            return IndependentCopula()
        else
            return new{typeof(θ)}(θ)
        end
    end
end

τ(C::MorgensternCopula) = 2*C.θ/9

function Distributions.cdf(C::MorgensternCopula, x::AbstractVector)
    θ = C.θ
    u1, u2 = x
    return u1*u2*(1+θ*(1-u1)*(1-u2))
end

function Distributions.pdf(C::MorgensternCopula, x::AbstractVector)
    θ = C.θ
    u1, u2 = x
    return 1+θ*(1-2*u1)*(1-2*u2)
end
function Distributions._rand!(rng::Distributions.AbstractRNG, C::MorgensternCopula, x::AbstractVector{T}) where {T<:Real}
    θ = C.θ
    u1, v2 = rand(rng, Distributions.Uniform(0,1),2)
    A = θ*(2*u1-1)-1
    B = (1-2*θ*(2*u1-1)+(θ^2)*(2*u1-1)^2 + 4*θ*v2*(2*u1-1))^(1/2)
    u2 = (2*v2)/(B-A)
    x[1] = u1
    x[2] = u2
    return x
end