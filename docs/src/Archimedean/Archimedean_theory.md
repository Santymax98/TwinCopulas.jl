# [Archimedean Copulas](@id Archimedean_theory)

*Archimedean copulas* are one of the most popular families within this field due to two fundamental characteristics: they have a closed and straightforward algebraic expression, and they exhibit a high level of symmetry. A copula $C_\varphi$, is an *Archimedean copula* if it has the functional form
$$C_\varphi(u_1, u_2) = \varphi(\varphi^{-1}(u_1) + \varphi^{-1}(u_2)),$$
where $\varphi$ is a continuous, strictly non-increasing function $\varphi: [0, \infty) \to [0, 1],$ with $\varphi(0)=1$ and $\lim_{x \to \infty}\varphi(x)=0,$ called the *Archimedean generator*.

In the context of bivariate Archimedean copulas, the generator $\varphi$ must also be completely monotone to ensure that $C_\varphi$ is a valid copula. A function $\varphi$ is said to be completely monotone if it satisfies $(-1)^k \varphi^{(k)}(x) \ge 0$ for all $x \ge 0$ and for all $k \ge 0$. This property guarantees that the corresponding copula exhibits appropriate dependency structures. To construct such generators, we can use the Williamson transform, which provides a method to generate completely monotone functions suitable for Archimedean copulas.

> **Definition (d-monotony [mcneil2009]):** A function $\varphi$ is said to be $d$-monotone if it has $d-2$ derivatives which satisfy 
>
> $(-1)^k \varphi^{(k)}(x) \ge 0 \;\forall k \in \{0, 1, \ldots, d-2\},$ and if $(-1)^{d-2}\varphi^{(d-2)}(x)$ is a non-increasing and convex function. 
>
> A function that is $d$-monotone for all $d$ is called **completely monotone**.

In the particular case when \(d = 2\), a function $\varphi$ is 2-monotone if it satisfies:

- $\varphi(x) \ge 0$
- $\varphi'(x) \le 0$ (i.e., $\varphi$ is non-increasing)
- $\varphi''(x) \ge 0$ (i.e., $\varphi'$ is non-decreasing and convex)

Therefore, a function that is 2-monotone is a non-increasing, convex function.

To summarize, a function $\varphi$ that is 2-monotone must be non-increasing and convex. If $\varphi$ is $d$-monotone for all $d$, it is called **completely monotone**.

In this package, there is an abstract type [`ArchimedeanCopula`] that provides a foundation for defining Archimedean copulas. Many Archimedean copulas are already implemented for you! See [the list of implemented Archimedean copulas] to get an overview.

If you do not find the one you need, you may define it yourself by subtyping `ArchimedeanCopula`. The API does not require much information, which is really convenient. Only the following methods are required:

* 𝘙(C::ArchimedeanCopula) =
* ϕ(C::ArchimedeanCopula, x) = 
* ϕ⁻¹(C::ArchimedeanCopula, x) =

That is, the Archimedean generator, its inverse, and its respective radial function. With the above, you can define a new copula in this package as follows: 

```julia
struct MyArchimedeanCopula{P} <: ArchimedeanCopula{P}
    θ::P
end

𝘙(C::MyArchimedeanCopula) = 1
ϕ(C::MyArchimedeanCopula, x) = exp(-x)
ϕ⁻¹(C::MyArchimedeanCopula, x) = -log(x)
```
It is not difficult to see that in this case $C_\varphi = \prod,$ showing that the independence copula $\prod$ is included in the Archimedean Family.

!!! info "Nomenclature information"
    We have called `𝘙()` the radial distribution, which is necessary for sampling a bivariate Archimedean copula. The letter we use is `\isansR`. In this document, we use $\varphi$, but in code, we use `\phi`.

# Advanced Concepts

Here, we present some important concepts of the theory of Archimedean copulas.

## Radial Distribution

The radial distribution, denoted as $𝘙(\cdot)$, is essential for generating samples from a bivariate Archimedean copula. It represents the distribution of the radial part of the copula and can be used to obtain the dependent structure between the two variables. The radial distribution is a key component in the sampling process and is used to generate the joint distribution from the marginal distributions.

For the class of Archimedean copulas, we propose three different algorithms to generate samples from a copula.

### Conditional Sampling for Bivariate Archimedean Copulas
To use this sampling method, we define the radial distribution as 1. This is represented in Julia as:

```Julia
𝘙(C::ArchimedeanCopula) = 1
```
When this is defined, TwinCopulas uses Algorithm 1, introduced in [reference], to sample from the copula and its derivative as follows. We adapt Algorithm 1 from [reference] for the case of bivariate Archimedean copulas. The input for the algorithm is a bivariate archimedean copula $C: [0,1]^2 \to [0,1].$

> **Algotithm 1: Bivariate Archimedean Copulas**
> 
> *(1)* Simulate $U_2 \sim \mathcal{U}[0,1]$
>
> *(2)* Compute (the right continuous version of) the function $$F_{U_1|U_2}(u_1)=\frac{\partial}{\partial u_2}C(u_1,u_2)=\varphi'(\varphi^{-1}(u_1)+\varphi^{-1}(u_2))\varphi^{-1'}(u_2), \quad u_2 \in [0,1].$$
> *(3)* Compute de generalized inverse of $F_{U_1|U_2},$ i.e $$F^{-1}_{U_1|U_2}(v)=\inf\{u_1 > 0: F_{U_1|U_2}(u_1)\geq v\}, \quad v \in [0,1].$$ 
> *(4)* Simulate $V \sim \mathcal{U}[0,1],$ independent of $U_2.$
>
> *(5)* Set $U_1 = F^{-1}_{U_1|U_2}(V)$ and return $(U_1, U_2).$

### Sampling from Archimedean Copula Using Frailty Distribution

### Algorithm

Here is a detailed algorithm for sampling from an Archimedean copula using the frailty distribution:

> **Algorithm 2: Sampling from Archimedean Copula Using Frailty Distribution**
> 
> *(1)* Sample a positive random variable $W$ with Laplace transform (Know Radial Distribution) $\varphi.$
>
> *(2)* Sample a random variable $Z \sim Erlang(2, 1)$ and compute $R=\frac{Z}{W}$ 
>
> *(3)* Sample i.i.d $E_1, E_2,$ where $E_k \sim Exp(1),$ and compute $S_k=\frac{E_k}{E_1+E_2}, \ \ k = 1,2.$ 
>
> *(4)* Return $(U_1, U_2)$ where $U_k=\varphi(RS_k), \ \ k=1,2.$

By using the frailty distribution approach, the algorithm ensures accurate and efficient sampling from the Archimedean copula, maintaining the desired dependence structure. For the copulas such that AMH, Clayton, Frank, Gumbel and Joe the literature shows its corresponding radial distribution.

```Julia
𝘙(C::GumbelCopula) = Stable(1 / C.θ)
```

### Sampling from Archimedean Copula Using Williamson Transform

In many cases the radial distribution is not kwon, in this case [`Mcneil`] solves this problem using the concept of d-monote generator together with the Williamson Transform.

**Definition (Williamson 2-Transform [mcneil2009]):** The *Williamson 2-Transform* of a positive random variable $R$ is defined as
$$\frak{W}_2F(R) = \mathbb{E}[\max(1-\frac{r}{R},0)], \ \ r \geq 0.$$

Furthermore, the distribution function $F_R$ of a postive random variable $R$ is completely characterized by Williamson transform of $\varphi$ archimedean generator by $\frak{W}_2^{-1}\varphi.$ That is
$$F_R(r)=1-\varphi(r)+r\varphi'(r).$$

This is sufficient to obtain the radial distribution necessary to obtain samples of a bivariate archimedean copula. This is very beneficial if you do not kwon the specific radial distribution for a bivariate archimedean copula. Follow the algorithm for this case:

> **Algorithm 3: Sampling from Archimedean Copula Using the Williamson 2-Transform**
> 
> *(1)* Sample a positive random variable $R\sim F(R)$ where $F(R)=\frak{W}_2^{-1}\varphi$.$
>
> *(2)* Sample i.i.d $E_1, E_2,$ where $E_k \sim Exp(1),$ and compute $S_k=\frac{E_k}{E_1+E_2}, \ \ k = 1,2.$ 
>
> *(3)* Return $(U_1, U_2)$ where $U_k=\varphi(RS_k), \ \ k=1,2.$

```Julia
𝘙(C::MyArchimedeanCopula) = RadialDist(C.θ)
```

We tested this implementation using `Nelsen2Copula` from TwinCopulas, this corresponds to copula number 2 from table 4.1 of Nelsen (2006). It works well for certain values of the parameter, for other cases there are numerical problems when obtaining the inverse distribution. We need more test to be able to generalize this function. For more details you can see [`Nelsen`]. 