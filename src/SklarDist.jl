struct SklarDist{C<:bicopula, M<:Tuple{<:Distributions.UnivariateDistribution, <:Distributions.UnivariateDistribution}} <: Distributions.ContinuousMultivariateDistribution
    copula::C
    margins::M

    function SklarDist(copula::C, margins::Vector{<:Distributions.UnivariateDistribution}) where {C<:bicopula}
        @assert length(margins) == 2 "Marginal distributions must be of length 2 for bivariate case"
        @assert all(margin -> margin isa Distributions.UnivariateDistribution, margins) "All margins must be univariate distributions"
        return new{C, Tuple{typeof(margins[1]), typeof(margins[2])}}(copula, (margins[1], margins[2]))
    end
end
# Definir la función length para SklarDist
Base.length(d::SklarDist) = 2

function Distributions.cdf(d::SklarDist, x::AbstractVector)
    u1 = Distributions.cdf(d.margins[1], x[1])
    u2 = Distributions.cdf(d.margins[2], x[2])
    return Distributions.cdf(d.copula, [u1, u2])
end

# Función para calcular la PDF conjunta
function Distributions.pdf(d::SklarDist, x::AbstractVector)
    u1 = Distributions.cdf(d.margins[1], x[1])
    u2 = Distributions.cdf(d.margins[2], x[2])
    du1 = Distributions.pdf(d.margins[1], x[1])
    du2 = Distributions.pdf(d.margins[2], x[2])
    copula_density = Distributions.pdf(d.copula, [u1, u2])
    return copula_density * du1 * du2
end

function Distributions._rand!(rng::Distributions.AbstractRNG, d::SklarDist, x::AbstractVector{T}) where {T<:Real}
    u1, u2 = Distributions.rand(rng, d.copula, 2)
    x[1] = Distributions.quantile(d.margins[1], u1)
    x[2] = Distributions.quantile(d.margins[2], u2)
    return x
end