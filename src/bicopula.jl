abstract type bicopula <: Distributions.ContinuousMultivariateDistribution end

# Función de longitud específica para bicopula
Base.length(::bicopula) = 2

function ρₛ(C::bicopula)
    F(x) = Distributions.cdf(C, x)
    z = zeros(2)
    i = ones(2)
    result, _ = Cubature.hcubature(F, z, i, reltol=sqrt(eps()))
    return 12 * result - 3
end

# Función para el tau de Kendall τ
function τ(C::bicopula)
    F(x) = Distributions.cdf(C, x)
    z = zeros(2)
    i = ones(2)
    result, _ = Cubature.hcubature(F, z, i, reltol=sqrt(eps()))
    return 4 * result - 1
end

function γ(C::bicopula)
    integrand1(u) = Distributions.cdf(C, [u, 1-u])
    integrand2(u) = u - Distributions.cdf(C, [u, u])
    
    integral1, _ = QuadGK.quadgk(integrand1, 0.0, 1.0)
    integral2, _ = QuadGK.quadgk(integrand2, 0.0, 1.0)
    
    return 4 * (integral1 - integral2)
end

# Función para calcular la medida de Blomqvist β
function β(C::bicopula)
    return 4 * Distributions.cdf(C, [0.5, 0.5]) - 1
end

# Función para calcular λ_u
function λᵤ(C::bicopula)
    u_vals = 0.999:0.00001:1.0
    λu_vals = [(Distributions.cdf(C, [u, u]) - 2u + 1) / (1 - u) for u in u_vals]
    return λu_vals[end]  # Aproximar el límite cuando u tiende a 1 desde la izquierda
end

# Función para calcular λ_l
function λₗ(C::bicopula)
    u_vals = 0.0:0.00001:0.001
    λl_vals = [Distributions.cdf(C, [u, u]) / u for u in u_vals if u != 0]
    return λl_vals[end]  # Aproximar el límite cuando u tiende a 0 desde la derecha
end

function dot_C(C::bicopula, x::AbstractVector)
    u1, u2 = x
    
    # Definimos una función que toma solo u2 como argumento
    cdf_wrt_u2(u2) = Distributions.cdf(C, [u1, u2])
    
    # Usamos ForwardDiff para calcular la derivada
    ForwardDiff.derivative(cdf_wrt_u2, u2)
end

function Distributions._rand!(rng::Distributions.AbstractRNG, C::bicopula, x::AbstractVector{T}) where {T<:Real}
    # Simular U2 y V uniformemente distribuidos
    u2, v = rand(rng, Distributions.Uniform(0, 1), 2)

    # Definir la función condicional
    function func(u)
        vectu = [u, u2]
        return dot_C(C, vectu) - v
    end
    
    u1 = Roots.find_zero(func, (eps(), 1 - eps()), Roots.Brent())
    
    # Asignar los valores calculados a x
    x[1] = u1
    x[2] = u2
    
    return x
end