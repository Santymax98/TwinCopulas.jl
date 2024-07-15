struct PStable{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    α::T
end

function PStable(α)
    if α <= 0 || α > 2
        throw(ArgumentError("Parameter α must be in (0, 2]"))
    end
    return PStable{typeof(α)}(α)
end
function Distributions.rand(rng::Distributions.AbstractRNG, d::PStable{T}) where T
    α = d.α
    if α == 1
        return 0.0
    else
        tα = 1 - α
        u = rand(rng, Distributions.Uniform(0,π))
        w = rand(rng, Distributions.Exponential())
        term1 = log((sin(tα*u)/w)^(tα/α))
        term2 = log(sin(α*u)/(sin(u)^(1/α)))
        return term1*term2
    end
end