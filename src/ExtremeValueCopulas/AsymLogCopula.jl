"""
    AsymLogCopula{P}

Fields:

  - α::Real - Dependency parameter
  - θ::Vector - Asymmetry parameters (size 2)

Constructor

    AsymLogCopula(α, θ)

The Asymmetric bivariate Logistic copula is parameterized by one dependence parameter ``\\alpha \\in [1, \\infty]`` and two asymmetry parameters ``\\theta_{i} \\in [0,1], i=1,2``. It is an Extreme value copula with Pickands dependence function: 

```math
A(t) = (\\theta_1^{\\alpha}(1-t)^{\\alpha} + \\theta_2^{\\alpha}t^{\\alpha})^{\\frac{1}{\\alpha}} + (\\theta_1 - \\theta_2)t + 1 - \\theta_1
```

References:
* Bivariate extreme value theory: models and estimation. Biometrika, 1988.
"""
struct AsymLogCopula{P} <: ExtremeValueCopula{P}
    α::P  # Parámetro de dependencia
    θ::Vector{P}  # Parámetros de asimetría (de tamaño 2 para el caso bivariado)
    function AsymLogCopula(α::P, θ::Vector{P}) where {P}
        if length(θ) != 2
            throw(ArgumentError("El vector θ debe tener 2 elementos para el caso bivariado"))
        elseif !(1 <= α)
            throw(ArgumentError("El parámetro α debe estar ser mayor o igual que 1"))
        elseif  !(0 <= θ[1] <= 1)  || !(0 <= θ[2] <= 1)  
            throw(ArgumentError("Todos los parámetros θ deben estar en el intervalo [0, 1]"))
        else
            return new{P}(α, θ)
        end
    end
end

# Definir la función A específica para la copula logística asimétrica bivariada
function 𝘈(C::AsymLogCopula, t::Real)
    α = C.α
    θ = C.θ
    
    A = ((θ[1]^α)*(1-t)^α + (θ[2]^α)*(t^α))^(1/α)+(θ[1]- θ[2])*t + 1 -θ[1]  
    return A
end