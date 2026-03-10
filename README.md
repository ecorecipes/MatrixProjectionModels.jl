# MatrixProjectionModels.jl

A Julia package for building and analyzing Matrix Projection Models (MPMs) — discrete-stage structured population models used in ecology and demography.

## Features

- **SciML-compatible**: `MPMProblem`/`solve()` pattern with `DirectIteration` and `EigenAnalysis` solvers
- **Full eigenanalysis**: asymptotic growth rate (λ), stable distribution, reproductive value, sensitivity, elasticity, damping ratio
- **Life history traits**: net reproductive rate, generation time, life expectancy, vital rates
- **Construction helpers**: Leslie and Lefkovitch matrices, mortality/fecundity models (Gompertz, Siler, Weibull, etc.), sparse transition syntax
- **Time-lagged models**: `LaggedMPM` with state augmentation (Kuss et al. 2008)
- **Stochastic models**: kernel-resampled and parameter-resampled stochasticity
- **COMPADRE integration**: load models from the COMPADRE Plant Matrix Database (via package extension)

## Quick Start

```julia
using MatrixProjectionModels

# 3-stage Leslie matrix
U = [0.0 0.0 0.0; 0.5 0.0 0.0; 0.0 0.3 0.0]
F = [0.0 3.0 1.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
mpm = MatrixProjectionModel(U, F)

println("λ = ", lambda(mpm))
println("Stable dist = ", stable_distribution(Matrix(mpm)))

# Population projection
n0 = [100.0, 50.0, 20.0]
prob = MPMProblem(mpm, n0, (0, 50))
sol = solve(prob, DirectIteration())
```

## Installation

This package is not yet registered in the Julia General registry. Install directly from GitHub (the [ProjectionModels.jl](https://github.com/ecorecipes/ProjectionModels.jl) dependency must be installed first):

```julia
using Pkg
Pkg.add(url="https://github.com/ecorecipes/ProjectionModels.jl")
Pkg.add(url="https://github.com/ecorecipes/MatrixProjectionModels.jl")
```

## Related

- [ProjectionModels.jl](https://github.com/ecorecipes/ProjectionModels.jl) — shared abstractions
- [IntegralProjectionModels.jl](https://github.com/ecorecipes/IntegralProjectionModels.jl) — continuous-state IPMs
- [CategoricalProjectionModels.jl](https://github.com/ecorecipes/CategoricalProjectionModels.jl) — categorical/functorial framework
- [COMPADRE](https://compadre-db.org/) — Plant Matrix Database
