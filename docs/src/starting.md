# Getting Started

## Installation

The *TwinCopulas* package is available through the Julia package system by running `Pkg.add("TwinCopulas")`.
Throughout, we assume that you have installed the package.

## Starting With a Bivariate Gaussian Copula

We start by drawing 1000 observations from a Gaussian Copula random variable.

The first step is to set up the environment:

```julia
julia> using Random, Distributions, TwinCopulas

julia> Random.seed!(123) # Setting the seed
```

Then, we create a bivariate Gaussian Copula distribution `G` and obtain samples using `rand`:

```julia
julia> G = GaussianCopula(0.5)
GaussianCopula{Float64}(θ=0.5)
```

The object `G` represents a probability distribution, in our case the bivariate Gaussian Copula.
One can query its properties such as the Kendall's tau ($\tau$, we need de letter \tau) :

```julia
julia> TwinCopulas.τ(G)
0.33333333333333337
```

We can also draw samples from `G` with `rand`.
```julia
julia> x = rand(G,1000)
2×1000 Matrix{Float64}:
 0.842684  0.151638  0.488056  0.0362257  0.474858  …  0.525896  0.640167  0.965114  0.0298152  0.632266
 0.195875  0.189835  0.777506  0.0273238  0.425615     0.41442   0.467543  0.533448  0.750643   0.359143
```

You can easily obtain the `pdf`, `cdf`, and many other functions for a distribution. For instance, the $\rho_s$ (Spearman's Rho), $\lambda_u, \lambda_l$ (Upper and lower tail dependence), $\beta$ (Blomqvist $\beta$) and $\gamma$ (Ginni's coefficient):

```julia
julia> cdf(G, [0.5,0.6])
0.3804363762931904

julia> pdf(G, [0.5,0.6])
1.1424140110638463
```

## Note: 

Functions such as `β` (Blomqvist's $\beta$), `τ` (Kendall's tau), `ρ_s` (Spearman's Rho), and tail dependence measures (`λ_u, λ_l`) must be prefixed with TwinCopulas because they are not directly exported by the module. For example, use TwinCopulas.β(copula) to access Blomqvist's $\beta$.

## Using Other Copulas

The package contains a large number of additional Copulas of two main types:

* `bicopula == ArrayLikeVariate{1}`
* `SklarDist == ArrayLikeVariate{1}`

Each type splits further into bivariate `Continuous`.

For instance, you can define the following Copulas (among many others):

```julia
julia> T = tCopula(ν, ρ) # bicopula (Elliptical)
julia> A = ClaytonCopula(θ)  # bicopula (Archimedean)
julia> E = GalambosCopula(θ) # bicopula (Extreme Value)
julia> R = ArchimaxCopula(A, E) # bicopula (Archimax)
```

In addition, you can create Copulas from univariate distributions (univariate margins):

```julia
julia> SklarDist(GaussianCopula(θ), [Normal(mu, sigma), Beta(α, β)])
```

To find out which parameters are appropriate for a given Copula `C`, you can use `fieldnames(D)`:

```julia
julia> fieldnames(tCopula)
(:θ, :ν)
julia> fieldnames(ArchimaxCopula)
(:Archimedean, :Extreme)
julia> fieldnames(SklarDist)
(:copula, :margins)
```

This tells you that a t Copula is initialized with dregree of freddom `ν` and  `θ`, ArchimaxCopula with an `Archimedean Copula` and `Extreme value Copula`, SklarDist need a bivariate `copula` (bicopula) and univariate distributions `margins`.