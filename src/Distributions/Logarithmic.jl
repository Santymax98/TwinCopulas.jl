struct Logarithmic{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    θ::T
    function Logarithmic(θ::Real) 
        if !(0 < θ < 1)
            throw(ArgumentError("Parameter θ must be in (0, 1)"))
        else
            return new{typeof(θ)}(θ)
        end 
    end
end

Statistics.mean(d::Logarithmic) = -d.θ / ((1 - d.θ) * log(1 - d.θ))
Statistics.var(d::Logarithmic) = -d.θ * (d.θ + log(1 - d.θ)) / ((1 - d.θ)^2 * (log(1 - d.θ))^2)

function Distributions.pdf(d::Logarithmic, x::Int)
    x > 0 ? -d.θ^x / (x * log(1 - d.θ)) : 0.0
end

function Distributions.cdf(d::Logarithmic, x::Int)
    θ = d.θ
    ks = 1:x
    s = sum(θ .^ ks ./ ks)
    return -s / log(1 - θ)
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::Logarithmic{T}) where T
    θ = d.θ
    h = log(1 - θ)
    u2 = rand(rng, Distributions.Uniform(0, 1))
    x = 1
    if u2 > θ
        return x
    else
        u1 = rand(rng, Distributions.Uniform(0, 1))
        q = 1 - exp(u1*h)
        if u2 < q^2
            return Int(trunc(1 + log(u2) / log(q)))
        else
            if u2 > q
                return 1
            else
                return 2
            end
        end
    end
end