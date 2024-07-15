struct AMHCopula{P} <: ArchimedeanCopula{P}
    θ::P
    function AMHCopula(θ)
        if !(-1 <= θ <= 1)
            throw(ArgumentError("El parámetro θ debe estar en [-1, 1]"))
        elseif θ == 0
            return IndependentCopula()
        elseif θ == 1
            return UtilCopula()
        else
            return new{typeof(θ)}(θ) 
        end
    end
end

𝘙(C::AMHCopula) = C.θ >= 0 ? 1 + Distributions.Geometric(1 - C.θ) : 1

ϕ(C::AMHCopula, x) = (1 - C.θ) / (exp(x) - C.θ)

ϕ⁻¹(C::AMHCopula, x) = log((1 - C.θ) / x + C.θ)

dϕ(C::AMHCopula, x) = -((1 - C.θ) * exp(x)) / (exp(x) - C.θ)^2

dϕ⁻¹(C::AMHCopula, x) = -(1 - C.θ) / (x * (C.θ * x + 1 - C.θ))

d²ϕ(C::AMHCopula, x) = ((1 - C.θ) * exp(x) * (exp(2 * x) - C.θ^2)) / (exp(x) - C.θ)^4

τ(C::AMHCopula) = (3 * C.θ - 2) / (3 * C.θ) - (2 * (1 - C.θ)^2 / (3 * C.θ^2)) * log(1 - C.θ)

ρₛ(C::AMHCopula) = (12 * (1 + C.θ) / C.θ^2) * PolyLog.reli2(1 - C.θ) - (24 * (1 - C.θ) / C.θ^2) * log(1 - C.θ) - (3 * (C.θ + 12) / C.θ)

λᵤ(C::AMHCopula) = 0.0

λₗ(C::AMHCopula) = 0.0