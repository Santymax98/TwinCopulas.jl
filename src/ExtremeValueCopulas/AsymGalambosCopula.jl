# Definición de la estructura AsymGalambosCopula
struct AsymGalambosCopula{P} <: ExtremeValueCopula{P}
    α::P  # Parámetro de dependencia
    θ::Vector{P}  # Parámetros de asimetría (de tamaño 2 para el caso bivariado)
    function AsymGalambosCopula(α::P, θ::Vector{P}) where {P}
        if length(θ) != 2
            throw(ArgumentError("El vector θ debe tener 2 elementos para el caso bivariado"))
        elseif !(0 <= α)
            throw(ArgumentError("El parámetro α debe estar ser mayor o igual que 0"))
        elseif  !(0 <= θ[1] <= 1)  || !(0 <= θ[2] <= 1)  
            throw(ArgumentError("Todos los parámetros θ deben estar en el intervalo [0, 1]"))
        elseif α == 0 || (θ[1] == 0 && θ[2] == 0)
            return IndependentCopula()
        elseif θ[1] == 1 && θ[2] == 1
            return GalambosCopula(α)
        else
            return new{P}(α, θ)
        end
    end
end

function 𝘈(C::AsymGalambosCopula, t::Real)
    α = C.α
    θ = C.θ

    term1 = (θ[1] * t)^(-α)
    term2 = (θ[2] * (1 - t))^(-α)
    
    inner_term = term1 + term2

    result = 1 - inner_term^(-1 / α)
    return result
end