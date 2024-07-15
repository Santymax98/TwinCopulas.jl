abstract type EllipticalCopula{P} <: bicopula end

function Distributions._rand!(rng::Distributions.AbstractRNG, C::EllipticalCopula, x::AbstractVector{T}) where {T<:Real}
    Σ = [1.0 C.θ; C.θ 1.0]
    A = LinearAlgebra.cholesky(Σ).L

    y = rand(rng, Distributions.Normal(), 2)
    if C isa GaussianCopula
    
        z = A*y
        x[1] = Distributions.cdf(Distributions.Normal(), z[1]) 
        x[2] = Distributions.cdf(Distributions.Normal(), z[2])

    elseif C isa tCopula

        w = rand(rng, Distributions.InverseGamma(C.ν/2.0, C.ν/2.0))
        z = sqrt(w)*A*y
        x[1] = Distributions.cdf(Distributions.TDist(C.ν), z[1]) 
        x[2] = Distributions.cdf(Distributions.TDist(C.ν), z[2])
    else
        throw(ArgumentError("Unsupported copula type"))
    end
    
    return x
end

τ(C::EllipticalCopula) = (2/π)*asin(C.θ)

ρₛ(C::EllipticalCopula) = (6/π)*asin(C.θ/2)

β(C::EllipticalCopula) = (2/π)*asin(C.θ)
