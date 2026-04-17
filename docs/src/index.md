# MatrixProjectionModels.jl

A Julia package for building and analyzing Matrix Projection Models (MPMs) -- discrete-stage structured population models used in ecology and demography.

## Features

- **SciML-compatible**: `MPMProblem`/`solve()` pattern with `DirectIteration` and `EigenAnalysis` solvers
- **Full eigenanalysis**: asymptotic growth rate (lambda), stable distribution, reproductive value, sensitivity, elasticity, damping ratio
- **Life history traits**: net reproductive rate, generation time, life expectancy, vital rates
- **Construction helpers**: Leslie and Lefkovitch matrices, mortality/fecundity models (Gompertz, Siler, Weibull, etc.)
- **Time-lagged models**: `LaggedMPM` with state augmentation
- **Stochastic models**: kernel-resampled and parameter-resampled stochasticity
- **COMPADRE integration**: load models from the COMPADRE Plant Matrix Database (via package extension)

## Quick Start

```julia
using MatrixProjectionModels

# 3-stage Leslie matrix
U = [0.0 0.0 0.0; 0.5 0.0 0.0; 0.0 0.3 0.0]
F = [0.0 3.0 1.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
mpm = MatrixProjectionModel(U, F)

println("lambda = ", lambda(mpm))
println("Stable dist = ", stable_distribution(Matrix(mpm)))

# Population projection
n0 = [100.0, 50.0, 20.0]
prob = MPMProblem(mpm, n0, (0, 50))
sol = solve(prob, DirectIteration())
```

## Installation

This package is not yet registered in the Julia General registry. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/ecorecipes/StructuredPopulationCore.jl")
Pkg.add(url="https://github.com/ecorecipes/MatrixProjectionModels.jl")
```

## Related Packages

- [StructuredPopulationCore.jl](https://github.com/ecorecipes/StructuredPopulationCore.jl) -- shared abstractions for projection models
- [IntegralProjectionModels.jl](https://github.com/ecorecipes/IntegralProjectionModels.jl) -- continuous-state integral projection models
- [CategoricalPopulationDynamics.jl](https://github.com/ecorecipes/CategoricalPopulationDynamics.jl) -- categorical/functorial framework
- [COMPADRE](https://compadre-db.org/) -- Plant Matrix Database
