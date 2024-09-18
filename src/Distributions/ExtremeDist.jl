struct ExtremeDist{C<:ExtremeValueCopula} <: Distributions.ContinuousUnivariateDistribution
    G::C
end

function Distributions.cdf(d::ExtremeDist, z)
    # Verificar que z esté en el soporte [0, 1]
    if z < 0 || z > 1
        return 0.0
    end
    copula = d.G
    return z + z * (1.0 - z) * (d𝘈(copula, z)/𝘈(copula, z))   # Se simplificó ya que la división se cancela
end

function _pdf(d::ExtremeDist, z)
    # Verificar que z esté en el soporte [0, 1]
    if z < 0 || z > 1
        return 0.0
    end
    
    copula = d.G
    A = 𝘈(copula, z)
    A_prime = d𝘈(copula, z)
    A_double_prime = d²𝘈(copula, z)
    return 1 + (1 - 2z) * A_prime / A + z * (1 - z) * (A_double_prime * A - A_prime^2) / A^2
end

function Distributions.quantile(d::ExtremeDist, p)
    if p < 0 || p > 1
        error("p debe estar en el rango [0, 1]")
    end

    cdf_func(x) = Distributions.cdf(d, x) - p
    copula = d.G

    # Automatically decide whether to use binary search or Brent
    if hasmethod(needs_binary_search, (typeof(copula),)) && needs_binary_search(copula)
        # Use binary search for copulas with large parameters     
        lower_bound = eps()
        upper_bound = 1.0 - eps()
        mid_point = (lower_bound + upper_bound) / 2
    
        while upper_bound - lower_bound > 1e-6  # Accuracy threshold
            mid_value = cdf_func(mid_point)
        
            if abs(mid_value) < 1e-6
                return mid_point
            elseif mid_value > 0
                upper_bound = mid_point
            else
                lower_bound = mid_point
            end
        
            mid_point = (lower_bound + upper_bound) / 2
        end
    
        return mid_point
    else
        # Use Brent for other copulations or if there are no problems
        return Roots.find_zero(cdf_func, (eps(), 1.0 - eps()), Roots.Brent())
    end
end

# Generar muestras aleatorias de la distribución radial usando la función cuantil
function Distributions.rand(rng::Distributions.AbstractRNG, d::ExtremeDist)
    u = rand(rng, Distributions.Uniform(0,1))  # Muestra de una distribución uniforme
    return Distributions.quantile(d, u)
end