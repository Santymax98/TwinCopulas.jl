"""
    MaresiasCopula{P}

Fields:

    - G::Function
    
Constructor

    MaresiasCopula(G)

The bivariate Maresias copula is parameterized by functions ``G, H: [0,1] \\to [0,1]``, such that ``G(u) = 2u - H(u), u \\in [0,1]``. It is constructed as: 

```math
C(u_1, u_2) = \\frac{1}{2}\\left ( G(u_1)G(u_2) + H(u_1)H(u_2) \\right )
```

References:
* [Mai2017](@cite) Simulating copulas: stochastic models, sampling algorithms, and applications. 2017.
"""
struct MaresiasCopula{G} <: bicopula
    G::G # debe ser una función
    function MaresiasCopula(G)
        if !is_valid_function(G)
            error("La función G debe ser creciente y estar entre 0 y 1")
        else
            return new{typeof(G)}(G)
        end
    end
end

function is_valid_function(G)
    u_values = range(0, stop=1, length=100)
    for u in u_values
        value = G(u)
        if value < 0 || value > 1
            return false
        end
    end
    return true
end

function Distributions.cdf(C::MaresiasCopula, x::AbstractVector)
    u1, u2 = x # solamente caso bivariado
    G = C.G
    H(u) = 2u - G(u) 
    prod_G = G(u1) * G(u2)
    prod_H = H(u1) * H(u2)
    return 0.5 * (prod_G + prod_H)
end

function Distributions.pdf(C::MaresiasCopula, x::AbstractVector)
    # Definimos una función que toma un vector y devuelve el cdf
    cdf_func = y -> Distributions.cdf(C, y)
    
    # Calculamos la derivada parcial mixta
    pdf_value = ForwardDiff.hessian(cdf_func, x)[1, 2]
    
    return pdf_value
end