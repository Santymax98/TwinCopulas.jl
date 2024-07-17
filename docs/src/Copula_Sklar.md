# Copulas

The law of 2-dimensional random vector $X = (X_1, X_2)$ defined on a probability space $(\Omega, \mathcal{F}, \mathbb{P}),$ is usually studied from its *distribution function*

$$F(x_1,x_2)= \mathbb{P}(X_1 \leq x_2, X_2 \leq x_2), \hspace{0.5cm} x_1, x_2 \in \mathbb{R}.$$

For $i = 1,2$ the distribution function $F_i$ of $X_i$ is called the (univariate) *marginal law* or *margin* and be retrieved from $F$ via

$$F_i(x_i) = \mathbb{P}(X_i \leq x_i)=F(x_i,\infty), \hspace{0.5cm} x_i \in \mathbb{R}.$$

It is important to note that knowing the marginal distributions $F_1, F_2$ is not sufficient to determine the joint distribution $F.$ Additionally, it is required to understand how the marginal distributions are coupled. This is achieved by means of a *copula* of $(X_1, X_2).$ Generally speaking, knowing the margins and a copula is equivalent to knowing the joint distribution. It is now appropriate to provide the definition of a copula.

> **Definition (Copula):**  
>
> *(1)* A function $C: [0, 1]^2 \to [0, 1]$ is called a *bivariate copula*, if there is a probability space $(\Omega, \mathcal{F}, \mathbb{P})$ supporting a random vector $(U_1, U_2)$ such that $U_k \sim \mathcal{U}[0,1]$ for all $k=1,2$ and
> 
> $$C(u_1,u_2)=\mathbb{P}(U_1 \leq u_1, U_2 \leq u_2), \hspace{0.5 cm} u_1, u_2 \in [0, 1].$$
>
> *(2)* On a probability space $(\Omega, \mathcal{F}, \mathbb{P})$ let $(U_1, U_2)$ be a random vector on $[0,1]^2$ whose joint distribution function is a copula
>
> $$C: [0,1]^2 \to [0,1].$$

```@example 1
using TwinCopulas
θ = 0.5 # Parameter
G = GaussianCopula(θ) # A 2-dimensional Gaussian Copula with parameter θ = 0.5.
```
This object is a random vector, and behaves exactly as you would expect a random vector from Distributions.jl to behave: you may sample it with rand(C,100), compute its pdf or cdf with pdf(C,x) and cdf(C,x), etc:

```@example 1
u = rand(G,10)
```
```@example 1
cdf(C,u)
```

> **Example (Independent Copula)**  
>
> The function $\Pi: [0,1]^2 \to [0,1],$ given by $$\Pi(u_1, u_2)=u_1u_2, \hspace{0.5 cm} u_1,u_2 \in [0,1],$$ is called the *Independence Copula.*
 
> **Example (Fréchet-Hoeffding bounds)**
>
>The Fréchet-Hoeffding bounds represent the extreme cases of dependence for a copula. These bounds are given by the following functions:
>
>- **Lower Fréchet-Hoeffding bound**:$$
  W(u_1, u_2) = \max\{u_1 + u_2 - 1, 0\}$$
>
>- **Upper Fréchet-Hoeffding bound**:$$M(u_1, u_2) = \min\{u_1, u_2\}$$
>
These limits can be viewed as the extreme cases of dependence in a copula:

- The **lower bound** $W(u_1, u_2)$ corresponds to the case of the most negative dependence, where the variables are perfectly anti-monotonic. This means that if one variable increases, the other decreases in a perfectly predictable manner.
  
- The **upper bound** $M(u_1, u_2)$ corresponds to the case of the most positive dependence, where the variables are perfectly monotonic. This means that if one variable increases, the other also increases in a perfectly predictable manner.

These bounds provide a way to understand the range of possible dependence structures that can be modeled by a copula, from complete negative dependence to complete positive dependence.

The syntax to use these copulas in TwinCopulas is the following:


```@example 2
using TwinCopulas
I = IndependentCopula() # A 2-dimensional Independent Copula.
M = MCopula() # A 2-dimensional upper bound.
W = WCopula() # A 2-dimensional lower bound.
```

# Sklar's Theorem
In probability theory and statistics, Sklar's Theorem is a fundamental result that connects multivariate distribution functions to their marginal distributions through a copula. This theorem provides a powerful framework for modeling and analyzing the dependence structure between random variables.

**Sklar's Theorem:**
Let $F$ a 2-dimensional distribution function with margins $F_1,F_2$. Then there exist a 2-dimensional copula $C$ such that for all $(x_1,x_2) \in \mathbb{R}^2$ it holds that 
$$F(x_1, x_2)=C(F_1(x_1),F_2(x_2)).$$  
If $F_1, F_2$ are continuous, then $C$ is unique. Conversely, if $C$ is a 2-dimensional copula and $F_1, F_2$ are univariate distribution functions, then he function $F$ is a 2-dimensional distribution function [`sklar1959`](@cite).



Sklar's Theorem is significant because it allows us to separate the marginal behavior of each variable from their dependence structure. The marginal distributions $F_1$ and $F_2$ describe the individual behavior of the random variables, while the copula $C$ captures how these variables are related to each other.

This separation is particularly useful in various applications:

1. **Modeling Flexibility**: By using copulas, we can model complex dependencies between variables while maintaining flexibility in choosing different marginal distributions. This is especially valuable in fields like finance, insurance, and risk management, where understanding dependencies between risks is crucial.

2. **Dependence Structure**: Copulas provide a way to quantify and visualize the dependence structure between variables, which is not possible by looking at marginals alone. This helps in understanding the nature and strength of the relationships between variables.

3. **Simulation and Sampling**: Sklar's Theorem facilitates the generation of multivariate data with specified marginals and dependence structures. This is useful for simulations and scenarios where realistic modeling of dependencies is required.

In summary, Sklar's Theorem is a cornerstone in the field of copula theory, enabling the decoupling of marginal distributions and dependence structures, and providing a versatile tool for statistical modeling and analysis of multivariate data.

The syntax for generating bivariate distributions with Sklar's theorem is as follows:

```@example 3
using TwinCopulas, Distributions
margins = [Normal(), Beta(3.5, 2.3)]
copula = ClaytonCopula(3.5)
F = SklarDist(copula, margins)
```

# Survival Copulas

In some cases, it is more convenient to describe the distribution of a random vector $(X_1, X_2)$ using its survival function instead of its distribution function, especially when the components $X_1$ and $X_2$ are interpreted as lifetimes. The survival function provides a more intuitive description in such contexts.

## Definition of the Bivariate Survival Function

For a bivariate random vector $(X_1, X_2)$ defined on a probability space $(\Omega, \mathcal{F}, \mathbb{P})$, the survival function $\overline{F}$ is defined as:

$$\overline{F}(x_1, x_2) := \mathbb{P}(X_1 > x_1, X_2 > x_2), \quad x_1, x_2 \in \mathbb{R}.$$

## Marginal Survival Functions

For $i = 1, 2$, the univariate marginal survival function $\overline{F}_i$ of $X_i$ can be retrieved from the bivariate survival function $\overline{F}$ as follows:

$$\overline{F}_i(x_i) = \mathbb{P}(X_i > x_i) = \overline{F}(x_i, -\infty) \quad \text{for} \; i = 1, \; \text{and} \; \overline{F}_i(x_i) = \overline{F}(-\infty, x_i) \quad \text{for} \; i = 2, \quad x_i \in \mathbb{R}.$$

## Relationship to Sklar's Theorem

Analogous to Sklar's Theorem for distribution functions, the survival function of a bivariate random vector can be decomposed into a copula and its marginal survival functions. This decomposition allows for a more nuanced understanding of the dependence structure between the variables.

## Survival Analog of Sklar's Theorem

Let $\overline{F}$ be a bivariate survival function with marginal survival functions $\overline{F}_1$ and $\overline{F}_2$. Then there exists a bivariate copula $\overline{C}$ such that for all $(x_1, x_2) \in \mathbb{R}^2$, the following holds:

$$\overline{F}(x_1, x_2) = \overline{C}(\overline{F}_1(x_1), \overline{F}_2(x_2)).$$

If the marginal survival functions $\overline{F}_1$ and $\overline{F}_2$ are continuous, then the copula $\overline{C}$ is unique. Conversely, if $\overline{C}$ is a bivariate copula and $\overline{F}_1$ and $\overline{F}_2$ are univariate survival functions, then the function $\overline{F}$ is a bivariate survival function.

The survival copula approach is particularly useful in reliability theory and survival analysis, where it is common to deal with the lifetimes of components or systems. The survival function provides a direct way to model the joint survival probabilities, and the survival copula captures the dependence structure between these lifetimes. This method allows for a clearer interpretation and more effective modeling of the joint behavior of the variables involved.

By utilizing survival copulas, one can effectively separate the marginal survival characteristics of each variable from their dependence structure, leading to a more flexible and comprehensive approach to multivariate survival analysis.

The syntax to use the survival copula in TwinCopulas is the following:


```@example 4
using TwinCopulas
cop = GumbelCopula(5.5)
S = SurvivalCopula(cop)
```

# Conditional Sampling for Bivariate Copulas

As the main sampling method, TwinCopulas implements conditional sampling based on a bivariate probability space, supporting two iid random variables. This algorithm is not restricted to any specific class of copulas, making it versatile for sampling arbitrary bivariate copulas. The only challenging step is calculating the partial derivative and its generalized inverse, which usually requires a well-defined analytical form of the copula and is not easy to obtain.

The imput for the algorithm is a bivariate copula $C: [0,1]^2 \to [0,1].$

> **Algotithm 1**
> 
> *(1)* Simulate $U_2 \sim \mathcal{U}[0,1]$
>
> *(2)* Compute (the right continuous version of) the function $$F_{U_1|U_2}(u_1)=\frac{\partial}{\partial u_2}C(u_1,u_2)|_{u_2=U_2}, \quad u_2 \in [0,1].$$
> *(3)* Compute de generalized inverse of $F_{U_1|U_2},$ i.e $$F^{-1}_{U_1|U_2}(v)=\inf\{u_1 > 0: F_{U_1|U_2}(u_1)\geq v\}, \quad v \in [0,1].$$ 
> *(4)* Simulate $V \sim \mathcal{U}[0,1],$ independent of $U_2.$
>
> *(5)* Set $U_1 = F^{-1}_{U_1|U_2}(V)$ and return $(U_1, U_2).$

