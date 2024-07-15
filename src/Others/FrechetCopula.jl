"""
    FrechetCopula{P}

Fields:

    - θ1::Real - parameter
    - θ2::Real - parameter
    
Constructor

    FrechetCopula(θ1, θ2)

The bivariate Fréchet copula is parameterized by ``\\theta_{i} \\in [0,1], i = 1,2``, such that ``\\theta_1 + \\theta_2 \\leq 1 \\``. It is constructed as: 

```math
C(u_1, u_2) = \\theta_1 \\min\\{u_1, u_2\\} + (1- \\theta_1 - \\theta_2)u_1u_2 + \\theta_2 \\max\\{u_1 + u_2 - 1, 0\\}
```

It has a few special cases: 
- When θ1 = θ2 = 0, it is the Independent Copula
- When θ1 = 1, it is the MCopula (Upper Frechet-Hoeffding bound) 
- When θ2 = 1, is is the WCopula (Lower Frechet-Hoeffding bound)

References:
* [nelsen2006](@cite) Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct FrechetCopula{P} <: bicopula
    θ1::P
    θ2::P

    function FrechetCopula(θ::Vararg{Real})
        if length(θ) !== 2
            throw(ArgumentError("FrechetCopula requires only 2 arguments."))
        end
        # Promover todos los parámetros al mismo tipo
        T = promote_type(typeof(θ[1]), typeof(θ[2]))
        θ1, θ2 = T(θ[1]), T(θ[2])
        
        # Validar que los parámetros sean no negativos
        if !(0 <= θ1 <= 1) || !(0 <= θ2 <= 1) 
            throw(ArgumentError("All θ parameters must be in [0,1]"))
        elseif !(θ1 + θ2 <= 1)
            throw(ArgumentError("θ1 + θ2 must be in [0,1]"))
        elseif θ1 == 1
            return MCopula()
        elseif θ2 == 1
            return WCopula()
        elseif θ1 == 0 && θ2 == 0
            return IndependentCopula()
        else
            return new{T}(θ1, θ2)
        end 
    end
end

function Distributions.cdf(C::FrechetCopula, u::AbstractVector)
    θ1, θ2 = C.θ1, C.θ2
    u1, u2 = u
    return θ1*min(u1, u2) + (1 - θ1 - θ2)*u1*u2 + θ2*max(u1 + u2 - 1, 0)
end

function Distributions.pdf(C::FrechetCopula, u::AbstractVector)
    throw(ArgumentError("There is no density for the Frechet Copula."))
end

function Distributions._rand!(rng::Distributions.AbstractRNG, C::FrechetCopula, x::AbstractVector{T}) where {T<:Real}
    p = [C.θ1, 1-C.θ1-C.θ2, C.θ2]
    
    u1, u2 = rand(rng, Distributions.Uniform(0,1), 2)
    
    z = rand(rng, Distributions.Categorical(p))
    
    if z == 1
        u1 = min(u1,u2)
        x[1] = u1
        x[2] = u1
    elseif z == 2
        x[1] = u1
        x[2] = u2
    else
        u1 = max(u1 + u2 - 1, 0)
        x[1] = u1
        x[2] = 1-u1
    end
    
    return x
end

function τ(C::FrechetCopula)
    θ1, θ2 = C.θ1, C.θ2
    term = (θ1 - θ2)(θ1 + θ2 + 2)
    return term/3
end