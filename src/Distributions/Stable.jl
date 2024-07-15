# Definición de la estructura Stable
struct Stable{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    α::T
end

# Constructor de la distribución estable
function Stable(α, β, γ, δ)
    if 0 < α <= 2 && -1 <= β <= 1 && γ > 0
        return Stable{typeof(α)}(α, β, γ, δ)
    else
        throw(ArgumentError("Parameters α must be in (0, 2], β must be in [-1, 1], and γ must be greater than 0"))
    end
end

function Distributions.rand(rng::Distributions.AbstractRNG, d::Stable{T}) where T
    θ = d.α
    U = rand(Distributions.Uniform(0,1))
    E = rand(Distributions.Exponential())
    U = π*(U - 0.5)
    term1 = sin(θ*(π/2 + U))
    term2 = cos(U)^(1/θ)
    term3 = (cos(U - θ*(π/2 + U))/E)^((1-θ)/θ)
    X = (term1/term2)*term3
    return X
end