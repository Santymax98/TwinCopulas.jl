"""
    UtilCopula{}

Constructor

    UtilCopula()

The bivariate UtilCopula is a simple copula that appears as a special case of several copulas, that has the form :  

```math
C(u_1, u_2) = \\frac{u_1u_2}{u_1+u_2 - u_1u_2}
```

It happends to be an Archimedean Copula, with generator : 

```math
\\phi(t) = 1 / (t + 1)
```

References:
* Nelsen, Roger B. An introduction to copulas. Springer, 2006.
"""
struct UtilCopula <: ArchimedeanCopula{Nothing} end

𝘙(C::UtilCopula) = 1

ϕ(C::UtilCopula, x) = 1 / (x + 1)

dϕ(C::UtilCopula, x) = -1 / (x + 1)^2

d²ϕ(C::UtilCopula, x) = 2 / (x + 1)^3

ϕ⁻¹(C::UtilCopula, x) = 1 / x - 1

dϕ⁻¹(C::UtilCopula, x) = -1 / x^2
