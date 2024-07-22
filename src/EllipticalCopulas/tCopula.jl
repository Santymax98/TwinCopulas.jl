"""
    tCopula{P}

Fields:
    - ν::Real - paremeter
    - θ::Real - Parameter 

Constructor

    tCopula(ν, θ)

The bivariate t copula. It is constructed as: 

```math
C(u_1, u_2; \\nu, \\theta) = t_{\\nu, \\theta}(t_{\\nu}^{-1}(u_1),t_{\\nu}^{-1}(u_2))
```
where ``t_{\\nu, \\theta}`` is the cumulative distribution function (CDF) of a bivariate t-distribution with ``\\nu \\in \\mathbb{R}^{+}`` degrees of freedom, zero means, and correlation ``\\theta \\in [-1, 1]``, and ``t_{\\nu}^{-1}`` is the quantile function of the standard t-distribution with ``\\nu`` degrees of freedom.

It has a few special cases:

- When θ = -1, it is the WCopula (Lower Frechet-Hoeffding bound)
- When θ = 0, it is the IndependentCopula
- When θ = 1, is is the MCopula (Upper Frechet-Hoeffding bound)

References:
* oe, Harry. Dependence modeling with Copulas. Chapman & Hall, 2014.
"""
struct tCopula{df, P} <: EllipticalCopula{P}
    θ::P   # Parámetro de correlación
    ν::df  # Grados de libertad

    function tCopula(ν::df, θ::P) where {df<:Real, P<:Real}
        if ν <= 0
            throw(ArgumentError("Los grados de libertad ν deben ser reales positivos"))
        end
        if !(-1 <= θ <= 1)
            throw(ArgumentError("El parámetro de correlación ρ debe estar entre -1 y 1"))
        elseif θ == 0
            return IndependentCopula()
        elseif θ == 1
            return MCopula()
        elseif θ == -1
            return WCopula()
        else
            return new{df, typeof(θ)}(θ, ν)
        end
    end
end

function Distributions.cdf(C::tCopula, x::AbstractVector)
    θ = C.θ
    ν = C.ν
    u1 = Distributions.quantile(Distributions.TDist(ν), x[1])
    u2 = Distributions.quantile(Distributions.TDist(ν), x[2])

    mvt_density = Distributions.MvTDist(ν, [0.0, 0.0], [1.0 θ; θ 1.0])

    integrand(v) = Distributions.pdf(mvt_density, v)
    lower_bounds = [-1000.0, -1000.0]
    upper_bounds = [u1, u2]
    
    result, _ = Cubature.hcubature(integrand, lower_bounds, upper_bounds)
    return result
end

function Distributions.pdf(C::tCopula, x::AbstractVector)
    θ = C.θ
    ν = C.ν
    u1 = Distributions.quantile(Distributions.TDist(ν), x[1])
    u2 = Distributions.quantile(Distributions.TDist(ν), x[2])
    num = (ν/2)*(SpecialFunctions.gamma(ν/2)^2)*(1 + (u1^2 + u2^2 -2*θ*u1*u2)/(ν*(1-θ^2)))^(-(ν+2)/2)
    dem = sqrt(1-θ^2)*(SpecialFunctions.gamma((ν+1)/2)^2)*((1+ (u1^2)/ν)*(1+(u2^2)/ν))^(-(ν+1)/2)
    return num/dem
end

function λᵤ(C::tCopula)
    θ = C.θ
    ν = C.ν
    term = -sqrt(((ν+1)*(1-θ))/(1+θ))
    return 2*Distributions.cdf(Distributions.TDist(ν+1), term)
end

function λₗ(C::tCopula)
    θ = C.θ
    ν = C.ν
    term = -sqrt(((ν+1)*(1-θ))/(1+θ))
    return 2*Distributions.cdf(Distributions.TDist(ν+1), term)
end