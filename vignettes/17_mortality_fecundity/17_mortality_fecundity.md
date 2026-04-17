# Mortality, Fecundity, and Error Estimation
Simon Frost

## Overview

This vignette covers three short topics:

1.  The abstract type hierarchies `AbstractMortalityModel` /
    `AbstractFecundityModel` for parametric vital-rate functions used in
    model construction.
2.  The `hx_to_px` / `px_to_hx` discrete hazard ↔ survival conversions
    underpinning life-table arithmetic.
3.  Bootstrap-based `calculate_errors` for matrix uncertainty.

## Setup

``` julia
using MatrixProjectionModels
using LinearAlgebra
using Random
```

## Vital-rate type hierarchy

``` julia
println("AbstractMortalityModel is abstract: ", isabstracttype(AbstractMortalityModel))
println("AbstractFecundityModel is abstract: ", isabstracttype(AbstractFecundityModel))
```

    AbstractMortalityModel is abstract: true
    AbstractFecundityModel is abstract: true

Concrete subtypes (e.g. `GompertzMortality`, `WeibullMortality`,
`LogisticFecundity`, `HadwigerFecundity`) all implement `model_survival`
or `model_fecundity` and the call-syntax `model(x)`.

``` julia
mort = GompertzMortality(0.001, 0.05)
fec  = LogisticFecundity(20.0, 5.0, 2.5)
println("typeof(mort) <: AbstractMortalityModel = ", typeof(mort) <: AbstractMortalityModel)
println("typeof(fec)  <: AbstractFecundityModel = ", typeof(fec)  <: AbstractFecundityModel)

ages = 0:5
println("hazards(0..5)   = ", round.(mort.(ages),  digits=4))
println("fecundity(0..5) = ", round.(fec.(ages),  digits=4))
```

    typeof(mort) <: AbstractMortalityModel = true
    typeof(fec)  <: AbstractFecundityModel = true
    hazards(0..5)   = [0.001, 0.0011, 0.0011, 0.0012, 0.0012, 0.0013]
    fecundity(0..5) = [0.0001, 0.0111, 1.5172, 18.4828, 19.9889, 19.9999]

## Discrete hazard ↔ survival

For discrete-time life tables the conversion between per-step survival
$p_x$ and per-step hazard $h_x = -\log p_x$ is linear-in-log:

``` julia
hx = [0.05, 0.07, 0.10, 0.20, 0.50]
px = hx_to_px(hx)
hx_back = px_to_hx(px)
println("px            = ", round.(px,      digits=4))
println("hx round-trip = ", round.(hx_back, digits=4))
println("max round-trip error = ", round(maximum(abs.(hx .- hx_back)), digits=8))
```

    px            = [0.9512, 0.9324, 0.9048, 0.8187, 0.6065]
    hx round-trip = [0.05, 0.07, 0.1, 0.2, 0.5]
    max round-trip error = 0.0

These compose with `lx_to_px`, `lx_to_hx`, `px_to_lx`, `hx_to_lx` to
convert any pair of life-table columns.

## Bootstrap error estimation

`calculate_errors(mpm, sample_size; type, n_boot)` resamples each entry
of the matrices under multinomial sampling error and reports either
standard errors (`:sem`) or 95% confidence intervals (`:ci95`).

``` julia
U = [0.4 0.0;
     0.5 0.6]
F = [0.0 2.0;
     0.0 0.0]
mpm = MatrixProjectionModel(U, F)

Random.seed!(42)
sem = calculate_errors(mpm, 100; type = :sem, n_boot = 200)
println("standard errors of A:")
for i in 1:size(sem.A, 1)
    println("  ", round.(sem.A[i, :], digits=4))
end
```

    standard errors of A:
      [0.0481, 0.1421]
      [0.0538, 0.0514]

``` julia
ci  = calculate_errors(mpm, 100; type = :ci95, n_boot = 200)
println("95% CI lower for A[2,2] = ", round(ci.lower[2, 2], digits=4))
println("95% CI upper for A[2,2] = ", round(ci.upper[2, 2], digits=4))
```

    95% CI lower for A[2,2] = 0.4995
    95% CI upper for A[2,2] = 0.69

## Summary

- `AbstractMortalityModel` and `AbstractFecundityModel` are the two
  parametric vital-rate hierarchies.
- `hx_to_px` and `px_to_hx` are the fundamental discrete hazard ↔
  survival conversions.
- `calculate_errors` is a bootstrap helper for matrix uncertainty
  reporting.
