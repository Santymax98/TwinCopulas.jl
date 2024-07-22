"""
    PlackettCopula{P}

Fields:
    - θ::Real - parameter

Constructor

    PlackettCopula(θ)

The bivariate Plackett Copula parameterized by ``\\theta > 0`` The [Plackett](https://www.cambridge.org/core/books/abs/copulas-and-their-applications-in-water-resources-engineering/plackett-copula/2D407DAB691623AB52CF74044B42C61F). It is constructed as:  

```math
C(u_1,u_2) = \\frac{1}{2}\\eta^{-1}(1+\\eta(u_1+u_2)-((1+\\eta(u_1+u_2))^2 - 4\\theta \\eta u_1 u_2)^{1/2}) 
```
where ``\\eta)=\\theta-1.``

It has a few special cases: 
- When θ = ∞, is is the MCopula (Upper Frechet-Hoeffding bound)
- When θ = 1, it is the IndependentCopula
- When θ = 0, is is the WCopula (Lower Frechet-Hoeffding bound) 

References:
* Joe, H. (2014). Dependence modeling with copulas. CRC press, Page.164
* Johnson, Mark E. Multivariate statistical simulation: A guide to selecting and generating continuous multivariate distributions. Vol. 192. John Wiley & Sons, 1987. Page 193.
* Nelsen, Roger B. An introduction to copulas. Springer, 2006. Exercise 3.38.
"""
struct PlackettCopula{P} <: bicopula
    θ::P  # Copula parameter

    function PlackettCopula(θ)
        if θ < 0
            throw(ArgumentError("Theta must be non-negative"))
        elseif θ == 0
            return MCopula()
        elseif θ == 1
            return IndependentCopula()
        elseif θ == Inf
            return WCopula()
        else
            return new{typeof(θ)}(θ)
        end
    end
end

# CDF calculation for bivariate Plackett Copula
function Distributions.cdf(C::PlackettCopula, x::AbstractVector)
    u1, u2 = x
    θ = C.θ
    η = θ - 1
    term1 = 1 + η * (u1 + u2)
    term2 = sqrt(term1^2 - 4 * θ * η * u1 * u2)
    return 0.5 * η^(-1) * (term1 - term2)
end

# PDF calculation for bivariate Plackett Copula
function Distributions.pdf(C::PlackettCopula, x::AbstractVector)
    u1, u2 = x
    θ = C.θ
    η = θ - 1
    term1 = θ * (1 + η * (u1 + u2 - 2 * u1 * u2))
    term2 = (1+η*(u1+u2))^2-4*(θ)*η*u1*u2
    return  (term1)*(term2)^(-3/2) 
end

function Distributions._rand!(rng::Distributions.AbstractRNG, C::PlackettCopula, x::AbstractVector{T}) where {T<:Real}
    θ = C.θ
    u1, t = rand(rng, Distributions.Uniform(0,1),2)
    a = t * (1 - t)
    b = θ + a * (θ - 1)^2
    cc = 2*a * (u1 * θ^2 + 1 - u1) + θ * (1 - 2a)
    d = sqrt(θ) * sqrt(θ + 4a * u1 * (1 - u1) * (1 - θ)^2)
    u2 = (cc - (1 - 2t) * d) / (2b)
    x[1] = u1
    x[2] = u2
    return x
end

# Calculate Spearman's rho based on the PlackettCopula parameters
function ρ(C::PlackettCopula)
    θ = C.θ
    return (θ+1)/(θ-1)-(2*θ*log(θ)/(θ-1)^2)
end