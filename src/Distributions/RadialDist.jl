struct RadialDist{C<:ArchimedeanCopula} <: Distributions.ContinuousUnivariateDistribution
    G::C
end

# Transformada de Williamson W_2
function 𝔚(G::ArchimedeanCopula, x)
    integrand(t) = (1 - x / t) * ϕ(G, t) * (t > x ? 1 : 0)
    integral, _ = QuadGK.quadgk(integrand, x, Inf)
    return integral
end

# Inversa de la transformada de Williamson para d = 2
function 𝔚⁻¹(G::ArchimedeanCopula, x)
    if  x < 0.0 
        return 0.0
    elseif x >= 1.0
        return 1.0
    else
        return 1 - ϕ(G, x) + x*dϕ(G, x)
    end
end

# CDF de la distribución radial
function Distributions.cdf(d::RadialDist, x)
    G = d.G
    if x <= 0
        return 0.0
    else
        return 𝔚⁻¹(G, x)
    end
end

# PDF de la distribución radial
function Distributions.pdf(d::RadialDist, x)
    G = d.G
    return (x * d²ϕ(G, x))
end

function Distributions.quantile(d::RadialDist, p)
    cdf_func(x) = Distributions.cdf(d, x) - p
    return Roots.find_zero(cdf_func, (eps(), 1000.0), Roots.Brent())
end

# Generar muestras aleatorias de la distribución radial usando la función cuantil
function Distributions.rand(rng::Distributions.AbstractRNG, d::RadialDist)
    u = rand(rng, Distributions.Uniform(0,1))  # Muestra de una distribución uniforme
    return Distributions.quantile(d, u)
end
