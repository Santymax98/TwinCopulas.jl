# TwinCopulas

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Santymax98.github.io/TwinCopulas.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Santymax98.github.io/TwinCopulas.jl/dev/)
[![Build Status](https://github.com/Santymax98/TwinCopulas.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Santymax98/TwinCopulas.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Build Status](https://api.cirrus-ci.com/github/Santymax98/TwinCopulas.jl.svg)](https://cirrus-ci.com/github/Santymax98/TwinCopulas.jl)
[![Coverage](https://codecov.io/gh/Santymax98/TwinCopulas.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Santymax98/TwinCopulas.jl)
[![Coverage](https://coveralls.io/repos/github/Santymax98/TwinCopulas.jl/badge.svg?branch=master)](https://coveralls.io/github/Santymax98/TwinCopulas.jl?branch=master)

*TwinCopulas* is a Julia package for probability distributions and associated functions, focusing on the implementation of copulas. This package includes:

* Probability density/mass functions (pdf).
* Sampling from a population or from a distribution.

## Installation

You can install *TwinCopulas* using Julia's package manager:

```julia
using Pkg
Pkg.add("TwinCopulas")
```

## Usage

Here is a basic example of how to use *TwinCopulas*:

```julia
Julia> using TwinCopulas
Julia> G = GaussianCopula(0.8)
```

## Roadmap

### Upcoming Features and Improvements

1. **Revised Calculation of Dependency Measures:**
   - Enhance the accuracy and efficiency of dependency measure calculations.

2. **Additional Copulas:**
   - Implement more copulas, potentially including dynamic copulas to capture time-varying dependencies.

3. **Copula Fitting (Parameter Estimation):**
   - Develop robust methods for parameter estimation of copulas to better fit data.

4. **Model Selection:**
   - Introduce model selection techniques to choose the best copula model for given data.

5. **Graphical Representations:**
   - Add functions to create informative plots and visualizations of copulas and their properties.

## Contributing

Contributions are welcome! Please open an issue if you have any questions, feature requests, or bug reports. Feel free to submit pull requests to improve the package.

## License

This project is licensed under the MIT License.

## Acknowledgements

We would like to thank all contributors and the Julia community for their support and contributions.
