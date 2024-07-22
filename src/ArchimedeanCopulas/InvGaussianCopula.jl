"""
    InvGaussianCopula{P}

Fields:

  - θ::Real - parameter

Constructor

    InvGaussianCopula(θ)

The bivariate Inverse Gaussian copula is parameterized by ``\\theta \\in [0,\\infty)``. It is an Archimedean copula with generator :

```math
\\phi(t) = \\exp{\\frac{1-\\sqrt{1+2θ^{2}t}}{θ}}.
```

More details about Inverse Gaussian Archimedean copula are found in :

    Mai, Jan-Frederik, and Matthias Scherer. Simulating copulas: stochastic models, sampling algorithms, and applications. Vol. 6. # N/A, 2017. Page 74.

It has a few special cases:
- When θ = 0, it is the IndependentCopula

References:
* Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct InvGaussianCopula{P} <: ArchimedeanCopula{P}
    θ::P
    function InvGaussianCopula(θ) 
        if !(0 <= θ) throw(ArgumentError("El parámetro θ debe estar en (0, ∞)")) 
        elseif θ == 0 
            return IndependentCopula()
        else 
            return new{typeof(θ)}(θ) 
        end 
    end
end

𝘙(C::InvGaussianCopula) = Distributions.InverseGaussian(C.θ, 1)

ϕ(C::InvGaussianCopula, x) = exp((1 - sqrt(1 + 2 * C.θ^2 * x)) / C.θ)

dϕ(C::InvGaussianCopula, x) = -C.θ / sqrt(1 + 2 * C.θ^2 * x) * exp((1 - sqrt(1 + 2 * C.θ^2 * x)) / C.θ)

d²ϕ(C::InvGaussianCopula, x) = C.θ^3 / (1 + 2 * C.θ^2 * x)^(3/2) * exp((1 - sqrt(1 + 2 * C.θ^2 * x)) / C.θ) + C.θ^2 / (1 + 2 * C.θ^2 * x) * exp((1 - sqrt(1 + 2 * C.θ^2 * x)) / C.θ)

ϕ⁻¹(C::InvGaussianCopula, x) = ((1 - C.θ * log(x))^2 - 1) / (2 * C.θ^2)

dϕ⁻¹(C::InvGaussianCopula, x) = - (1 - C.θ * log(x)) / (x * C.θ)

λᵤ(C::InvGaussianCopula) = 0.0

λₗ(C::InvGaussianCopula) = 0.0
