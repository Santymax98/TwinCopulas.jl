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

function Distributions.rand(rng::Distributions.AbstractRNG, d::Sibuya{T}) where {T <: Real}
    u = rand(rng, T)
    if u <= d.θ
        return T(1)
    end
    xMax = 1/eps(T)
    Ginv = ((1-u)*SpecialFunctions.gamma(1-d.θ))^(-1/d.θ)
    fGinv = floor(Ginv)
    if Ginv > xMax 
        return fGinv
    end
    if 1-u < 1/(fGinv*SpecialFunctions.beta(fGinv,1-d.θ))
        return ceil(Ginv)
    end
    return fGinv
end