"""
    EmpiricalCopula{M} 

Fields:

    -data::M - 2xn matrix observations 

Constructor

    EmpiricalCopula(data)

Let ``(x_k, y_k), k = 1,2, \\ldots, n`` denote a sample of size ``n`` from continuous bivariate distribution. The Empirical Copula is the function

```math
C_n(\\frac{i}{n},\\frac{j}{n}) = \\frac{\\text{number of pairs (x,y) in the sample with} x \\leq x_{(i)}, y \\leq y_{(j)}}{n}
```
For more details see

References:
* Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct EmpiricalCopula{M} <: bicopula
    data::M
    pseudos::M

    function EmpiricalCopula(data)
        @assert size(data, 1) == 2 "Data must have two dimensions for bivariate case"
        pseudo_obs = pseudos(data)
        new{typeof(data)}(data, pseudo_obs)
    end
end

# Función para calcular la CDF empírica
function Distributions.cdf(C::EmpiricalCopula, x::AbstractVector)
    u1, u2 = clamp.(x, 0, 1)
    n = size(C.data, 2)
    count = sum((C.pseudos[1, :] .<= u1) .& (C.pseudos[2, :] .<= u2))
    return count / n
end

function Distributions.pdf(C::EmpiricalCopula, x::AbstractVector)
    n = size(C.data, 2)
    u1, u2 = clamp.(x, 0, 1)
    h = 1 / sqrt(n)  # Ancho de banda adaptativo
    kernel_vals = epanechnikov_kernel.((C.pseudos[1, :] .- u1) / h) .* epanechnikov_kernel.((C.pseudos[2, :] .- u2) / h)
    return sum(kernel_vals) / (n * h^2)
end


function frequency(C::EmpiricalCopula, x::AbstractVector)
    n = size(C.data, 2)
    i = round(Int, x[1] * n)
    j = round(Int, x[2] * n)
    return sum((C.pseudos[1, :] .== i/n) .& (C.pseudos[2, :] .== j/n)) / n
end
# Función para generar muestras
function Distributions._rand!(rng::Distributions.AbstractRNG, C::EmpiricalCopula, x::AbstractVector{T}) where {T<:Real}
    n = size(C.data, 2)
    
    # Generar U(0,1) aleatorios
    u1, u2 = rand(rng), rand(rng)
    
    # Encontrar los índices más cercanos
    i1 = searchsortedfirst(C.pseudos[1, :], u1)
    i2 = searchsortedfirst(C.pseudos[2, :], u2)
    
    # Interpolación lineal simple para u1
    if i1 == 1
        x[1] = C.pseudos[1, 1]
    elseif i1 > n
        x[1] = C.pseudos[1, n]
    else
        x[1] = C.pseudos[1, i1-1] + (u1 - C.pseudos[1, i1-1]) * (C.pseudos[1, i1] - C.pseudos[1, i1-1]) / (C.pseudos[1, i1] - C.pseudos[1, i1-1])
    end
    
    # Interpolación lineal simple para u2
    if i2 == 1
        x[2] = C.pseudos[2, 1]
    elseif i2 > n
        x[2] = C.pseudos[2, n]
    else
        x[2] = C.pseudos[2, i2-1] + (u2 - C.pseudos[2, i2-1]) * (C.pseudos[2, i2] - C.pseudos[2, i2-1]) / (C.pseudos[2, i2] - C.pseudos[2, i2-1])
    end
    
    return x
end

function pseudos(sample::AbstractArray)
    d, n = size(sample)
    @assert d == 2 "Sample must have exactly two dimensions for bivariate case"
    ranks = zeros(eltype(sample), d, n)
    for i in 1:d
        ranks[i, :] = (StatsBase.ordinalrank(sample[i, :]) .- 0.5) ./ n
    end
    return ranks
end

function epanechnikov_kernel(u)
    return abs(u) <= 1 ? 0.75 * (1 - u^2) : 0.0
end