```@meta
CurrentModule = TwinCopulas
```

# TwinCopulas

The [*TwinCopulas*](https://github.com/Santymax98/TwinCopulas.jl) package provides an extensive collection of bivariate copulas and related tools for their manipulation and analysis. In particular, *TwinCopulas* implements:

* **Archimedean Copulas**: such as Clayton's Copula, Gumbel's Copula, Frank's Copula, and more.
* **Elliptical Copulas**: such as the Gaussian Copula and t-Student Copula.
* **Extreme Value Copulas**: such as the Galambos Copula, Hüsler-Reiss Copula, and others.
* **Sklar Distribution**: Implements Sklar's Theorem, which states that any bivariate joint distribution can be expressed in terms of its marginal distributions and a copula that describes the dependency structure between the variables.
* **Other Copulas**: including Archimax Copulas, Empirical Copulas, and more.

* **Tools for Copula Analysis**: functions for sample generation, evaluation of the cumulative distribution function (CDF) and the probability density function (PDF), calculation of dependency measures, and more.

This package is designed to facilitate working with bivariate copulas in Julia, providing an intuitive and efficient interface for teaching and understanding important concepts in this area. If you need another package to investigate and work with high-dimensional data, you can see and use [*Copulas*](https://github.com/lrnv/Copulas.jl) by Oskar Laverny and Santiago Jiménez (me).

*TwinCopulas* is specifically designed for academic environments because working in two dimensions is simpler and clearer than working in higher dimensions. This is particularly beneficial in statistical inference when maximizing functions to obtain maximum likelihood estimators, or in mathematical statistics when seeking closed-form solutions to equations. Additionally, the package aims to implement a wide range of important concepts in copula theory and dependency analysis, most of which are well-documented in the literature for two dimensions. However, this package is also suitable for generating new knowledge and intensive work, as it is natively implemented in Julia. By following the conventions of [*Distributions*](https://github.com/JuliaStats/Distributions.jl), *TwinCopulas* is compatible with many other packages such as [*Turing*](https://github.com/TuringLang/Turing.jl), [*StatsBase*](https://github.com/JuliaStats/StatsBase.jl), and others, making it a versatile tool for both teaching and research.
