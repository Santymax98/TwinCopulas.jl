"""
    AsymMixedCopula{P}

Fields:

  - θ::Vector - parameters (size 2)

Constructor

    AsymMixedCopula(θ)

The Asymmetric bivariate Mixed copula is parameterized by two parameters ``\\theta_{i}, i=1,2`` that must meet the following conditions:
* θ₁ ≥ 0
* θ₁+θ₂ ≤ 1
* θ₁+2θ₂ ≤ 1
* θ₁+3θ₂ ≥ 0

It is an Extreme value copula with Pickands dependence function: 

```math
A(t) = \\theta_{2}t^3 + \\theta_{1}t^2-(\\theta_1+\\theta_2)t+1
```

It has a few special cases:

- When θ₁ = θ₂ = 0, it is the Independent Copula
- When θ₂ = 0, it is the Mixed Copula

References:
* [Tawn1988](@cite) Bivariate extreme value theory: models and estimation. Biometrika, 1988.
"""
struct AsymMixedCopula{P} <: ExtremeValueCopula{P}
    θ::Vector{P}  # Parámetros de asimetría (de tamaño 2 para el caso bivariado)

    function AsymMixedCopula(θ::Vector{P}) where {P}
        if length(θ) != 2
            throw(ArgumentError("El vector θ debe tener 2 elementos para el caso bivariado"))
        elseif !(0 <= θ[1])
            throw(ArgumentError("El parámetro θ₁ debe estar ser mayor o igual que 0"))
        elseif  !(θ[1]+θ[2] <= 1) 
            throw(ArgumentError("la suma de θ₁+θ₂ ≤ 1"))
        elseif !(θ[1]+2*θ[2] <= 1)
            throw(ArgumentError("la suma de θ₁+2θ₂ ≤ 1"))
        elseif !(0 <= θ[1]+3*θ[2])
            throw(ArgumentError("la suma de 0 ≤ θ₁+3θ₂"))
        elseif θ[1] == 0 && θ[2] == 0
            return IndependentCopula()
        else
            return new{P}(θ)
        end
    end
end

# Definir la función A específica para la copula logística asimétrica bivariada
function 𝘈(C::AsymMixedCopula, t::Real)
    θ = C.θ
    
    A = θ[2]*t^3 + θ[1]*t^2-(θ[1]+θ[2])*t+1
    return A
end