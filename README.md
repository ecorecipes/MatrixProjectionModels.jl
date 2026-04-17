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

## Vignettes

| # | Vignette | Description |
|---|----------|-------------|
| 1 | [Introduction to Matrix Projection Models](https://github.com/ecorecipes/MatrixProjectionModels.jl/blob/main/vignettes/01_introduction/01_introduction.md) | Core concepts: constructing MPMs, eigenanalysis, population projection |
| 2 | [Age-Structured (Leslie) Models](https://github.com/ecorecipes/MatrixProjectionModels.jl/blob/main/vignettes/02_leslie_models/02_leslie_models.md) | Building Leslie matrices from age-specific survival and fecundity |
| 3 | [Model Construction and Random Generation](https://github.com/ecorecipes/MatrixProjectionModels.jl/blob/main/vignettes/03_construction/03_construction.md) | Constructing MPMs from mortality and fecundity models (Gompertz, Siler, etc.) |
| 4 | [Vital Rate Extraction and Decomposition](https://github.com/ecorecipes/MatrixProjectionModels.jl/blob/main/vignettes/04_vital_rates/04_vital_rates.md) | Extracting survival, growth, and reproduction rates from projection matrices |
| 5 | [Life Tables and Age-from-Stage Analysis](https://github.com/ecorecipes/MatrixProjectionModels.jl/blob/main/vignettes/05_life_tables/05_life_tables.md) | Constructing life tables and age-from-stage distributions |
| 6 | [Life History Traits](https://github.com/ecorecipes/MatrixProjectionModels.jl/blob/main/vignettes/06_life_history/06_life_history.md) | Net reproductive rate, generation time, life expectancy, longevity |
| 7 | [Perturbation Analysis](https://github.com/ecorecipes/MatrixProjectionModels.jl/blob/main/vignettes/07_perturbation/07_perturbation.md) | Sensitivity, elasticity, and LTRE decomposition |
| 8 | [Stochastic Matrix Projection Models](https://github.com/ecorecipes/MatrixProjectionModels.jl/blob/main/vignettes/08_stochastic/08_stochastic.md) | Kernel-resampled and parameter-resampled stochastic models |
| 9 | [Population Simulation and Density Dependence](https://github.com/ecorecipes/MatrixProjectionModels.jl/blob/main/vignettes/09_simulation/09_simulation.md) | Individual-based simulation and density-dependent models |
| 10 | [Comparative Demography with COMPADRE and COMADRE](https://github.com/ecorecipes/MatrixProjectionModels.jl/blob/main/vignettes/10_comparative/10_comparative.md) | Loading and analyzing models from the COMPADRE/COMADRE databases |
| 11 | [Sparse Transition Constructors](https://github.com/ecorecipes/MatrixProjectionModels.jl/blob/main/vignettes/11_transitions/11_transitions.md) | Building matrices from sparse transition specifications |
| 12 | [Time-Lagged Matrix Projection Models](https://github.com/ecorecipes/MatrixProjectionModels.jl/blob/main/vignettes/12_time_lag/12_time_lag.md) | State-augmented models with delayed fecundity (Kuss et al. 2008) |

## Installation

This package is not yet registered in the Julia General registry. Install directly from GitHub (the [StructuredPopulationCore.jl](https://github.com/ecorecipes/StructuredPopulationCore.jl) dependency must be installed first):

```julia
using Pkg
Pkg.add(url="https://github.com/ecorecipes/StructuredPopulationCore.jl")
Pkg.add(url="https://github.com/ecorecipes/MatrixProjectionModels.jl")
```

## Related

- [StructuredPopulationCore.jl](https://github.com/ecorecipes/StructuredPopulationCore.jl) — shared abstractions
- [IntegralProjectionModels.jl](https://github.com/ecorecipes/IntegralProjectionModels.jl) — continuous-state IPMs
- [FiniteStatePopulationDynamics.jl](https://github.com/ecorecipes/FiniteStatePopulationDynamics.jl) — discrete-state continuous-time dynamics
- [ContinuousStatePopulationDynamics.jl](https://github.com/ecorecipes/ContinuousStatePopulationDynamics.jl) — continuous-state continuous-time dynamics
- [CategoricalPopulationDynamics.jl](https://github.com/ecorecipes/CategoricalPopulationDynamics.jl) — categorical/functorial framework
- [PhysiologicallyBasedDemographicModels.jl](https://github.com/ecorecipes/PhysiologicallyBasedDemographicModels.jl) — application-level PBDM reference suite
- [COMPADRE](https://compadre-db.org/) — Plant Matrix Database
