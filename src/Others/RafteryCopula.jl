"""
    RafteryCopula{P}

Fields:

    - θ::Real - parameter
    
Constructor

    RafteryCopula(θ)

The bivariate Raftery Copula is parameterized by ``\\theta \\in [0,1]``. It is constructed as: 

```math
C(u_1, u_2) = \\min\\{u_1,u_2\\} +\\frac{1-\\theta}{1+\\theta}(u_1u_2)^{1/(1-\\theta)}\\left \\{1-(\\max\\{u_1,u_2\\})^{-(1+\\theta)/(1-\\theta)}  \\right \\}
```
It has a few special cases: 
- When θ = 0, it is the Independent Copula
- When θ = 1, it is the MCopula (Upper Frechet-Hoeffding bound)

References:
* Joe, Harry, Multivariate Models and Multivariate Dependence Concepts, Chapman & Hall. 1997.
"""
struct RafteryCopula{P} <: bicopula
    θ::P  # Copula parameter

    function RafteryCopula(θ)
        if !(0 <= θ <= 1)
            throw(ArgumentError("Theta must be in [0,1]"))
        elseif θ == 0
            return IndependentCopula()
        elseif θ == 1
            return MCopula()
        else
            return new{typeof(θ)}(θ)
        end
    end
end

# CDF calculation for bivariate Plackett Copula
function Distributions.cdf(C::RafteryCopula, x::AbstractVector)
    u1, u2 = x
    θ = C.θ
    term1 = ((1-θ)/(1+θ))*(u1*u2)^(1/(1-θ))
    term2 = 1-(max(u1,u2))^(-(1+θ)/(1-θ))
    return min(u1,u2)+(term1*term2)
end

# PDF calculation for bivariate Plackett Copula
function Distributions.pdf(C::RafteryCopula, x::AbstractVector)
    θ = C.θ
    u1, u2 = sort(x)
    term1 = 1/(θ^2 - 1)
    term2 = -1-θ*(u2)^((-θ-1)/(1-θ))
    term3 = (u1*u2)^(θ/(1-θ))
    return term1*term2*term3
end

function Distributions._rand!(rng::Distributions.AbstractRNG, C::RafteryCopula, x::AbstractVector{T}) where {T<:Real}
    θ = C.θ
    u, u1, u2 = rand(rng, Distributions.Uniform(0,1),3)
    j = rand(rng, Distributions.Bernoulli(θ))
    
    x[1] = (u1^(1 - θ)) * (u)^j
    x[2] = (u2^(1 - θ)) * (u)^j
    return x
end

# CCorregir este rho
function ρ(C::RafteryCopula)
    θ = C.θ
    return (θ+1)/(θ-1)-(2*θ*log(θ)/(θ-1)^2)
end