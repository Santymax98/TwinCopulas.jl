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
#aux functions
function G(d::Sibuya, x::Real)
    θ = d.θ
    x_min = (gamma(1 - θ))^(-1/θ)
    if x < x_min
        throw(ArgumentError("x must be greater than or equal to $x_min."))
    end
    return 1 - (1 / (x^θ * gamma(1 - θ)))
end

function G⁻¹(d::Sibuya, y::Real)
    if y < 0 || y > 1
        throw(ArgumentError("y must be in [0, 1]."))
    end
    θ = d.θ
    return ((1 - y) * gamma(1 - θ))^(-1/θ)
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::Sibuya{T}) where T
    θ = d.θ
    u = rand(rng, Uniform(0, 1))
    if u <= θ
        return 1
    else
        G_inv_u = G⁻¹(d, u)
        if Distributions.cdf(d, floor(G_inv_u)) <= u
            return ceil(G_inv_u)
        else
            return floor(G_inv_u)
        end
    end
end