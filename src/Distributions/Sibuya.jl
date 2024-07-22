struct Sibuya{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    θ::T
    function Sibuya(θ::Real)
        if !(0 < θ <=1)
            throw(ArgumentError("Parameter θ must be in (0, 1]"))
        else
            new{typeof(θ)}(θ)
        end
    end
end

#Statistics.mean(d::Sibuya) = 0.0
#Statistics.var(d::Sibuya) = 0.0

function Distributions.pdf(d::Sibuya, x::Int)
    x > 0 ? (-1)^(x+1)*(1/(x*SpecialFunctions.beta(x, d.θ-x+1))) : 0.0
end

function Distributions.cdf(d::Sibuya, x::Int)
    θ = d.θ
    if x < 0
        throw(ArgumentError("x must be a positive integer"))
    else
        return 1-(1/(x*SpecialFunctions.beta(x,1-θ)))
    end
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::Sibuya{T}) where T
    θ = d.θ
    Z1 = rand(rng, Distributions.Exponential())
    Z2 = rand(rng, Distributions.Gamma(1 - θ))
    Z3 = rand(rng, Distributions.Gamma(θ))
    
    λ = (Z1 * Z2) / Z3
    
    # Limitar λ utilizando eps
    xMax = 1 / eps(T)
    if λ > xMax
        λ = xMax
    end

    # Verificar el tipo de dato antes de pasar a Poisson
    if λ > typemax(Int32)
        return 1 + rand(rng, Distributions.Poisson(trunc(Int64, λ)))
    else
        return 1 + rand(rng, Distributions.Poisson(trunc(Int32, λ)))
    end
end