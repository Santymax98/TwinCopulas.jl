struct ArchimaxCopula{A<:ArchimedeanCopula, E<:ExtremeValueCopula} <: bicopula
    Archimedean::A
    Extreme::E
    function ArchimaxCopula(Archimedean::A, Extreme::E) where {A<:ArchimedeanCopula, E<:ExtremeValueCopula}
        new{A, E}(Archimedean, Extreme)
    end
end


function Distributions.cdf(C::ArchimaxCopula, x::AbstractVector)
    Archimedean = C.Archimedean
    Extreme = C.Extreme
    term1 = ϕ⁻¹(Archimedean, x[1]) + ϕ⁻¹(Archimedean, x[2])
    term2 = (ϕ⁻¹(Archimedean, x[2]))/term1
    term3 = 𝘈(Extreme, term2)
    return ϕ(Archimedean, term1*term3)
end

function Distributions.pdf(C::ArchimaxCopula, x::AbstractVector)
    # Definimos una función que toma un vector y devuelve el cdf
    cdf_func = y -> Distributions.cdf(C, y)
    
    # Calculamos la derivada parcial mixta
    pdf_value = ForwardDiff.hessian(cdf_func, x)[1, 2]
    
    return pdf_value
end

function τ(C::ArchimaxCopula)
    Archimedean = C.Archimedean
    Extreme = C.Extreme
    tA = τ(Extreme)
    tphi = τ(Archimedean)
    return tA + (1 - tA)*tphi
end

function Distributions._rand!(rng::Distributions.AbstractRNG, C::ArchimaxCopula, x::AbstractVector{T}) where {T<:Real}
    Archimedean = C.Archimedean
    Extreme = C.Extreme
    
    v1, v2 = Distributions.rand(rng, Extreme, 2)
    M = rand(rng, 𝘙(Archimedean))
    
    x[1] = ϕ(Archimedean, -log(v1)/M)
    x[2] = ϕ(Archimedean, -log(v2)/M)
    
    return x
end