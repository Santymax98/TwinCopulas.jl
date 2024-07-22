"""
    GaussianCopula{P}

Fields:

  - θ::Real - Parameter 

Constructor

    GaussianCopula(θ)

The bivariate [Gaussian](https://en.wikipedia.org/wiki/Copula_(probability_theory)#Gaussian_copula) copula. It is constructed as: 

```math
C(u_1, u_2; \\theta) = \\Phi_{\\theta}(\\Phi^{-1}(u_1),\\Phi^{-1}(u_2))
```
where ``\\Phi_{\\theta}`` is the cumulative distribution function (CDF) of a standard bivariate normal distribution with correlation coefficient ``\\theta \\in 
[-1, 1]`` and ``\\Phi^{-1}``is the quantile function of the standard normal distribution.

It has a few special cases:

- When θ = -1, it is the WCopula (Lower Frechet-Hoeffding bound)
- When θ = 0, it is the IndependentCopula
- When θ = 1, is is the MCopula (Upper Frechet-Hoeffding bound)

References:
* Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct GaussianCopula{P} <: EllipticalCopula{P}
    θ::P # Copula parameter
    function GaussianCopula(θ)
        if !(-1 <= θ <= 1)
            throw(ArgumentError("Theta must be in (-1, 1)"))
        elseif θ == 0
            return IndependentCopula()
        elseif θ == -1
            return WCopula()
        elseif θ == 1
            return MCopula()
        else
            return new{typeof(θ)}(θ)
        end
    end
end

function Distributions.cdf(C::GaussianCopula, x::AbstractVector)
    θ = C.θ
    u1 = Distributions.quantile(Distributions.Normal(), x[1])
    u2 = Distributions.quantile(Distributions.Normal(), x[2])

    mvn_density = Distributions.MvNormal([0.0, 0.0], [1.0 θ; θ 1.0])

    integrand(v) = Distributions.pdf(mvn_density, v)
    lower_bounds = [-1000.0, -1000.0]
    upper_bounds = [u1, u2]
    
    result, _ = Cubature.hcubature(integrand, lower_bounds, upper_bounds)
    return result
end

function Distributions.pdf(C::GaussianCopula, x::AbstractVector)
    θ = C.θ
    u1 = Distributions.quantile(Distributions.Normal(), x[1])
    u2 = Distributions.quantile(Distributions.Normal(), x[2])
    factor1 = 1 / sqrt(1 - θ^2)
    exponent = (2θ * u1 * u2 - θ^2 * (u1^2 + u2^2)) / (2 * (1 - θ^2))
    density = factor1 * exp(exponent)
    return density
end