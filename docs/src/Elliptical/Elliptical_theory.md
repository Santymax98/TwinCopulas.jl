# [Elliptical Copulas](@id Elliptical_theory)

An *Elliptical Copula* is defined as the copula of the related elliptical distribution $F$. Its analytical form is obtained via *Sklar's theorem* from the distribution function $F$. A copula $C$ is an *Elliptical Copula* if it has the functional form
$$C(u_1, u_2) = F(F_1^{-1}(u_1), F_2^{-1}(u_2)), \quad u_1, u_2 \in [0, 1]^2,$$
where $F_k^{-1}$ are the univariate quantile functions, $k=1,2$.

*Elliptical copulas* are widely used and well-documented in the literature due to their unique characteristics, such as their ability to model symmetric dependencies and their flexibility in capturing various dependency structures. Unlike extreme value copulas and Archimedean copulas, elliptical copulas do not have a closed-form expression, making them distinct in their formulation and application.

Two well-known examples of elliptical copulas are the Gaussian copula and the t copula. These copulas are particularly important due to their widespread use in finance and risk management:

> **Bivariate Gaussian Copula**: The bivariate normal or Gaussian copula $C_{gaussian}$ is the copula of $(X, Y) \sim \mathcal{N}_2(0, P),$ where $P$ is a correlation matrix. The functional form is obtained by $$C_{Gaussian}(u_1,u_2)=\Phi_2(\Phi^{-1}(u_1),\Phi^{-1}(u_2)), \quad u_1,u_2 \in [0,1]^2,$$ where $\Phi_2$ is the joint distribution function of $(X,Y),$ and $\Phi^{-1}$ is the quantile function of the univariate standard normal distribution.

> **Bivariate t Copula**: The bivariate $t$-Copula $C_{t,\nu}$ is the copula of $(X,Y)\sim t_2(0,P,\nu),$ where $P$ is a correlation matrix. The analytical form is obtained by $$C_{t,\nu}(u_1,u_2)=t_{2,\nu}(t_{\nu}^{-1}(u_1),t_{\nu}^{-1}(u_2)), \quad u_1,u_2 \in [0,1]^2,$$ where $t_{2,\nu}$ is the joint distribution function of $(X,Y)$ and $t_{\nu}^{-1}$ is the quantile function of the univariate standar $t$-distribution with $\nu$ degrees of freedom.

Elliptical copulas, including the Gaussian and t copulas, are powerful tools for modeling dependencies between variables, particularly when the dependency structure is symmetric and not overly influenced by tail behavior.

In this package, there is an abstract type [`EllipticalCopula`](@ref) that provides a foundation for defining elliptical copulas. The Gaussian and t copulas are already implemented for you! You can utils these as follow

```Julia
julia> G = GaussianCopula(0.8)
GaussianCopula{Float64}(θ=0.8) #Gaussian Copula with param 0.8

julia> T = tCopula(2, -0.3)
tCopula{Int64, Float64}(θ=-0.3, ν=2) #t Copula with 2 degree of freedom and param -0.3
```

You may want to use an elliptical copula customized for your needs. TwinCopulas allows the implementation of your own elliptic copula. You can implement it in the following way

```julia
struct MyEllipticalCopula{P} <: EllipticalCopula{P}
    θ::P #param o params need for created a correlation matrix
end
```

In this way you create a structure of your copula and you can implement the `CDF`, `PDF` or whatever functions you want following the conventions used in this package

```Julia
function Distributions.cdf(C::MyEllipticalCopula, x::AbstractVector)
    cdf_value = #Implement the specific cdf of your copula
    return cdf_value
end

function Distributions.pdf(C::MyElliptical, x::AbstractVector)
    pdf_value = #Implement the specific pdf of your copula
    return pdf_value
end
```

## Sampling for Bivariate Elliptical Copulas

In the same way as in other classes of domes, to generate samples of a Gaussian copula we can use the conditional method, you only need the first derivative of the copula function, however in many cases this can be complex. You could use ForwardDiff for an approximation.

In *TwinCopulas* we use specific methods for each elliptical copula. Below we detail these algorithms.

### Simulating the Gaussian Copula

Imput the bivariate Gaussian Copula with param $\theta$.

> *(1)* We form the correlation matrix $\Sigma = \begin{pmatrix}
1.0 & \theta\\ 
\theta & 1.0
\end{pmatrix}$ for the bivariate case.
>
> *(2)* Compute the Cholesky descomposition of $\Sigma,$ providing a lower triangular matrix $A,$ satisfying $AA'=\Sigma$
>
> *(3)* Simulate $Y_1, Y_2 \sim \mathcal{N}(0,1)$
>
> *(4)* Compute the random vector $$\begin{pmatrix}
X_1\\ 
X_2
\end{pmatrix}=A\cdot \begin{pmatrix}
Y_1\\ 
Y_2
\end{pmatrix}$$
>which has bivariate normal distribution with zero mean vector, unit variances, and correlation matrix $\Sigma$
>
> *(5)* Return $(U_1, U_2),$ where $U_i=\Phi(X_i), i=1,2,$ where  $\Phi$ is the CDF of standard normal distribution.

### Simulating the Gaussian Copula

Imput the bivariate $t$-Copula with params $(\nu,\theta)$.

> *(1)* We form the correlation matrix $\Sigma = \begin{pmatrix}
1.0 & \theta\\ 
\theta & 1.0
\end{pmatrix}$ for the bivariate case.
>
> *(2)* Compute the Cholesky descomposition of $\Sigma,$ providing a lower triangular matrix $A,$ satisfying $AA'=\Sigma$
>
> *(3)* Simulate $Y_1, Y_2 \sim \mathcal{N}(0,1)$ and simulate a rando variable $W\sim inv\Gamma(\frac{\nu}{2}, \frac{\nu}{2})$
>
> *(4)* Compute the random vector $$\begin{pmatrix}
X_1\\ 
X_2
\end{pmatrix}=\sqrt{W}A\cdot \begin{pmatrix}
Y_1\\ 
Y_2
\end{pmatrix}$$
>
> *(5)* Return $(U_1, U_2),$ where $U_i=t_{\nu}(X_i), i=1,2,$ where  $t_{\nu}$ is the CDF of standard $t$-distribution.