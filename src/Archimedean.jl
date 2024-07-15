abstract type ArchimedeanCopula{P} <: bicopula end

𝘙(C::ArchimedeanCopula) = 𝘙(C.θ)

ϕ(C::ArchimedeanCopula, x) = throw(ArgumentError("Function ϕ must be defined for specific copula"))

dϕ(C::ArchimedeanCopula, x) = throw(ArgumentError("Function dϕ must be defined for specific copula"))

d²ϕ(C::ArchimedeanCopula, x) = throw(ArgumentError("Function d²ϕ must be defined for specific copula"))

ϕ⁻¹(C::ArchimedeanCopula, x) = throw(ArgumentError("Function ϕ⁻¹ must be defined for specific copula"))

dϕ⁻¹(C::ArchimedeanCopula, x) = throw(ArgumentError("Function dϕ⁻¹ must be defined for specific copula"))

function Distributions.cdf(C::ArchimedeanCopula, u::Vector)
    inner_term = ϕ⁻¹(C, u[1]) + ϕ⁻¹(C, u[2])
    return ϕ(C, inner_term)
end

function Distributions.pdf(C::ArchimedeanCopula, u::Vector)
    u1, u2 = u
    inner_term = ϕ⁻¹(C, u1) + ϕ⁻¹(C, u2)
    product = dϕ⁻¹(C, u1) * dϕ⁻¹(C, u2)
    return d²ϕ(C, inner_term) * product
end

function Distributions._rand!(rng::Distributions.AbstractRNG, C::ArchimedeanCopula, x::AbstractVector{T}) where {T<:Real}
    if 𝘙(C) == 1
        u2, v = rand(rng, Distributions.Uniform(0, 1), 2)
        func(u) = (u - ϕ⁻¹(C, u)/dϕ⁻¹(C, u)) - v
    
        if func(eps()) > 0.0
            u1 = 0.0
        else
            u1 = Roots.find_zero(func, (eps(), 1-eps()), Roots.Brent())
        end    
        x[1] = ϕ(C, u2 * ϕ⁻¹(C, u1))
        x[2] = ϕ(C, (1 - u2) * ϕ⁻¹(C, u1))
        return x
    elseif 𝘙(C) == RadialDist(C)
        Y = rand(rng, Distributions.Exponential(), 2)
        S = Y/sum(Y)
        R = rand(rng, 𝘙(C))
        x[1] = ϕ(C, R*S[1])
        x[2] = ϕ(C, R*S[2])
        return x
    else
        W = Distributions.rand(rng, 𝘙(C))
        Z = rand(rng, Distributions.Erlang(2))
        R = Z / W
        E = rand(rng, Distributions.Exponential(), 2)
        S = sum(E)
        
        x[1] = ϕ(C, R * E[1] / S)
        x[2] = ϕ(C, R * E[2] / S)
        return x
    end
end

function τ(C::ArchimedeanCopula)
    integrand(x) = ϕ⁻¹(C, x) / dϕ⁻¹(C, x)
    result, _ = QuadGK.quadgk(integrand, 0, 1)
    return 1 + 4 * result
end
