"""
    MixedCopula{P}

Fields:

    - θ::Real - parameter
    
Constructor

    MixedCopula(θ)

The bivariate Mixed copula is parameterized by ``\\alpha \\in [0,1]``. It is an Extreme value copula with Pickands dependence function: 

```math
A(t) = \\theta t^2 - \\theta t + 1
```

It has a few special cases: 
- When θ = 0, it is the IndependentCopula

References:
* Bivariate extreme value theory: models and estimation. Biometrika, 1988.
"""
struct MixedCopula{P} <: ExtremeValueCopula{P}
    θ::P  # Parámetro de la copula
    function MixedCopula(θ)
        if !(0 <= θ <= 1)
            throw(ArgumentError("El parámetro θ debe estar en el intervalo [0, 1]"))
        elseif θ == 0
            return IndependentCopula()
        else
            return new{typeof(θ)}(θ)
        end
    end
end

𝘈(C::MixedCopula, t::Real) = C.θ*t^2 - C.θ*t + 1