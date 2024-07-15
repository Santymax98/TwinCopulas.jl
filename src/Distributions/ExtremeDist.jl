struct ExtremeDist{C<:ExtremeValueCopula} <: Distributions.ContinuousUnivariateDistribution
    G::C
end

function Distributions.cdf(d::ExtremeDist, z)
    copula = d.G
    return z + z*(1 - z)*(d𝘈(copula, z)/𝘈(copula, z)) 
end
using Distributions
function _pdf(d::ExtremeDist, z)
    copula = d.G
    A = 𝘈(copula, z)
    A_prime = d𝘈(copula, z)
    A_double_prime = d²𝘈(copula, z)
    return 1 + (1 - 2z) * A_prime / A + z * (1 - z) * (A_double_prime * A - A_prime^2) / A^2
end

function Distributions.quantile(d::ExtremeDist, p)
    cdf_func(x) = Distributions.cdf(d, x) - p
    return Roots.find_zero(cdf_func, (eps(), 1-eps()), Roots.Brent())
end

# Generar muestras aleatorias de la distribución radial usando la función cuantil
function Distributions.rand(rng::Distributions.AbstractRNG, d::ExtremeDist)
    u = rand(rng, Distributions.Uniform(0,1))  # Muestra de una distribución uniforme
    return Distributions.quantile(d, u)
end