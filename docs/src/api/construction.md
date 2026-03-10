# Construction

Tools for constructing matrix projection models from parametric mortality and fecundity models, Leslie/Lefkovitch builders, and sampling utilities.

## Mortality Models

Parametric mortality models used to generate age-specific survival schedules.

```@docs
AbstractMortalityModel
GompertzMortality
GompertzMakehamMortality
ExponentialMortality
SilerMortality
WeibullMortality
WeibullMakehamMortality
model_survival
```

## Fecundity Models

Parametric fecundity models used to generate age-specific fertility schedules.

```@docs
AbstractFecundityModel
LogisticFecundity
StepFecundity
VonBertalanffyFecundity
NormalFecundity
HadwigerFecundity
model_fecundity
```

## Leslie & Lefkovitch Builders

Convenience functions for constructing and randomly generating Leslie and Lefkovitch matrices.

```@docs
make_leslie_mpm
rand_lefko_mpm
rand_leslie_set
rand_lefko_set
```

## Sampling & Error

Functions for adding sampling error to projection matrices and computing confidence intervals.

```@docs
add_mpm_error
calculate_errors
compute_ci
```
