# Life Table Response Experiments (LTRE)

## Overview

A **Life Table Response Experiment (LTRE)** decomposes differences or variation in the population growth rate $\lambda$ into contributions from individual matrix elements (vital rates). While sensitivity and elasticity analysis asks "what *would* happen if a vital rate changed?", LTRE analysis asks "what *did* contribute to the observed difference?"

This vignette demonstrates four LTRE approaches available in MatrixProjectionModels.jl:

1. **Classical fixed-design LTRE** (Caswell 2001) — sensitivity-weighted element differences
2. **Classical random-design LTRE** — covariance-based decomposition of Var(λ)
3. **Exact fixed-design LTRE** (Hernandez et al. 2023) — fANOVA decomposition with no linear approximation
4. **Exact random-design LTRE** — exact decomposition of Var(λ) into main effects + interactions
5. **Small-noise approximation LTRE** (Davison et al. 2019) — decomposes stochastic growth rate differences

## Setup

```@example mpm
using MatrixProjectionModels
using LinearAlgebra
using Plots
```

## Example: Loggerhead Sea Turtle

We use a simplified 5-stage life cycle for loggerhead sea turtles (*Caretta caretta*), comparing a "baseline" population with a population where headstarting increases juvenile survival.

```@example mpm
# Baseline matrix (Crouse et al. 1987, simplified)
A_baseline = [
    0.0   0.0   0.0   0.0   61.9
    0.675 0.703 0.0   0.0   0.0
    0.0   0.047 0.657 0.0   0.0
    0.0   0.0   0.019 0.682 0.0
    0.0   0.0   0.0   0.061 0.809
]

# Treatment: headstarting increases egg-to-juvenile survival
A_headstart = [
    0.0   0.0   0.0   0.0   61.9
    0.85  0.703 0.0   0.0   0.0   # increased from 0.675 to 0.85
    0.0   0.047 0.657 0.0   0.0
    0.0   0.0   0.019 0.750 0.0   # increased from 0.682 to 0.75
    0.0   0.0   0.0   0.061 0.809
]

println("λ baseline   = ", round(lambda(A_baseline), digits=4))
println("λ headstart  = ", round(lambda(A_headstart), digits=4))
println("Δλ           = ", round(lambda(A_headstart) - lambda(A_baseline), digits=4))
```

## Classical Fixed-Design LTRE

The classical approach approximates:
$$\lambda_{\text{treat}} - \lambda_{\text{ref}} \approx \sum_{i,j} (a_{ij}^{\text{treat}} - a_{ij}^{\text{ref}}) \times s_{ij}$$

where $s_{ij}$ is the sensitivity evaluated at the midpoint matrix.

```@example mpm
result = ltre(A_headstart, A_baseline)

println("Sum of contributions = ", round(sum(result.contributions), digits=5))
println("Actual Δλ            = ", round(result.delta_lambda, digits=5))
println("Approximation error  = ", round(abs(sum(result.contributions) - result.delta_lambda), digits=6))
```

The sum of contributions closely approximates Δλ. Let's visualize which transitions drive the difference:

```@example mpm
stages = ["Egg", "Sm Juv", "Lg Juv", "Subadult", "Adult"]
C = result.contributions

# Find non-zero contributions
significant = findall(abs.(C) .> 1e-6)
for idx in significant
    i, j = Tuple(idx)
    println("  $(stages[j]) → $(stages[i]): contribution = $(round(C[i,j], digits=5))")
end
```

```@example mpm
heatmap(C, 
    xticks=(1:5, stages), yticks=(1:5, stages),
    xlabel="From stage", ylabel="To stage",
    title="LTRE Contributions to Δλ",
    color=:RdBu, clim=(-0.05, 0.05))
```

## Exact Fixed-Design LTRE

The exact method (Hernandez et al. 2023) uses functional ANOVA to decompose Δλ **without** the linear (sensitivity) approximation. It captures interaction effects between vital rates.

```@example mpm
exact_result = exact_ltre(A_headstart, A_baseline)

println("Number of varying elements: ", length(exact_result.indices_varying))
println("Number of effects (main + interactions): ", length(exact_result.effects))
println("Sum of ALL effects = ", round(sum(exact_result.effects), digits=10))
println("Actual Δλ          = ", round(lambda(A_headstart) - lambda(A_baseline), digits=10))
```

The exact method recovers Δλ *perfectly* — no approximation error. Let's examine main effects vs interactions:

```@example mpm
# Main effects are those with single indices
main_effects = [(exact_result.effect_indices[i], exact_result.effects[i]) 
                for i in 1:length(exact_result.effects)
                if length(exact_result.effect_indices[i]) == 1]

interactions = [(exact_result.effect_indices[i], exact_result.effects[i]) 
                for i in 1:length(exact_result.effects)
                if length(exact_result.effect_indices[i]) > 1]

println("Main effects:")
for (idx, eff) in main_effects
    eff == 0 && continue
    println("  Element $(idx[1]): $(round(eff, digits=6))")
end

println("\nInteraction effects:")
for (idx, eff) in interactions
    abs(eff) < 1e-10 && continue
    println("  Elements $idx: $(round(eff, digits=6))")
end

println("\nSum main effects:  ", round(sum(e for (_, e) in main_effects), digits=6))
println("Sum interactions:  ", round(sum(e for (_, e) in interactions), digits=6))
```

## Classical Random-Design LTRE

When we have multiple annual matrices (not a treatment vs. reference), a random-design LTRE decomposes Var(λ) into contributions from variance and covariance in matrix elements.

```@example mpm
# Simulate 5 years of environmental variation
A_years = [
    A_baseline .* (1 .+ 0.05 .* randn(5, 5)),
    A_baseline .* (1 .+ 0.05 .* randn(5, 5)),
    A_baseline .* (1 .+ 0.05 .* randn(5, 5)),
    A_baseline .* (1 .+ 0.05 .* randn(5, 5)),
    A_baseline .* (1 .+ 0.05 .* randn(5, 5)),
]
# Ensure non-negative
A_years = [max.(A, 0) for A in A_years]

rresult = ltre_random(A_years)

println("Observed Var(λ) = ", round(rresult.var_lambda, digits=6))
println("Sum of contributions = ", round(sum(rresult.contributions), digits=6))
println("λ of mean matrix = ", round(rresult.lambda_mean, digits=4))
```

The diagonal of the contribution matrix gives variance contributions; off-diagonal gives covariance contributions:

```@example mpm
# Top variance contributors (diagonal)
diag_contribs = diag(rresult.contributions)
n2 = length(diag_contribs)
n = Int(sqrt(n2))
top_idx = sortperm(abs.(diag_contribs), rev=true)[1:5]
println("Top 5 variance contributors:")
for idx in top_idx
    j = div(idx - 1, n) + 1
    i = mod(idx - 1, n) + 1
    println("  a[$i,$j] ($(stages[j])→$(stages[i])): $(round(diag_contribs[idx], digits=8))")
end
```

## Exact Random-Design LTRE

The exact random-design LTRE provides fANOVA-based decomposition of Var(λ):

```@example mpm
er_random = exact_ltre(A_years)

println("Number of varying elements: ", length(er_random.indices_varying))
println("Number of effects: ", length(er_random.effects))
println("Sum of effects = ", round(sum(er_random.effects), digits=8))
```

## Small-Noise Approximation LTRE (SNA-LTRE)

When comparing two sets of annual matrices, the SNA-LTRE (Davison et al. 2019) decomposes the difference in stochastic growth rate into four components:

1. **Mean**: shifts in mean vital rates
2. **Elasticity**: shifts in elasticity structure
3. **CV**: shifts in temporal variation
4. **Correlation**: shifts in temporal correlations

```@example mpm
# "Good" environment: lower variation
A_good = [A_baseline .* (1 .+ 0.02 .* randn(5, 5)) for _ in 1:8]
A_good = [max.(A, 0) for A in A_good]

# "Bad" environment: higher variation
A_bad = [A_baseline .* (1 .+ 0.10 .* randn(5, 5)) .* 0.9 for _ in 1:8]
A_bad = [max.(A, 0) for A in A_bad]

sna = sna_ltre(A_good, A_bad)

println("SNA-LTRE decomposition:")
println("  r_good (log λ_s) = ", round(sna.r_treatment, digits=4))
println("  r_bad  (log λ_s) = ", round(sna.r_reference, digits=4))
println("  Δr               = ", round(sna.delta_r, digits=4))
println()
println("  Component contributions:")
println("    Mean:        ", round(sum(sna.cont_mean), digits=5))
println("    Elasticity:  ", round(sum(sna.cont_elasticity), digits=5))
println("    CV:          ", round(sum(sna.cont_cv), digits=5))
println("    Correlation: ", round(sum(sna.cont_correlation), digits=5))
```

## Comparison: Classical vs Exact

A key advantage of the exact LTRE is that it has **zero approximation error** and accounts for interactions:

```@example mpm
# Classical
classical = ltre(A_headstart, A_baseline)
classical_error = abs(sum(classical.contributions) - classical.delta_lambda)

# Exact
exact = exact_ltre(A_headstart, A_baseline)
exact_error = abs(sum(exact.effects) - (lambda(A_headstart) - lambda(A_baseline)))

println("Classical LTRE approximation error: ", round(classical_error, digits=10))
println("Exact LTRE approximation error:     ", round(exact_error, digits=10))
println()
println("The exact method captures interaction effects that the classical")
println("method misses, which is especially important when multiple vital")
println("rates change simultaneously.")
```

## The G-Matrix

The exact LTRE uses the G-matrix (Poelwijk et al. 2016) to convert "responses" into "effects" via Möbius inversion:

```@example mpm
G = gmatrix(3)
println("G-matrix for 3 varying parameters (8×8):")
display(G)
```

## Summary

| Method | Decomposes | Approximation | Interactions |
|--------|-----------|---------------|-------------|
| Classical fixed | Δλ | First-order (sensitivity) | No |
| Exact fixed | Δλ | None (exact) | Yes |
| Classical random | Var(λ) | First-order | Pairwise only |
| Exact random | Var(λ) | None (exact) | All orders |
| SNA-LTRE | Δlog(λ_s) | Second-order | Mean/elas/CV/corr |
| Stochastic LTRE | Δlog(λ_s) | Simulation-based | Mean + SD |

## References

- Caswell, H. (2001) *Matrix Population Models*, Ch. 10. Sinauer.
- Hernandez, C.M. et al. (2023) An exact method for LTRE. *Methods Ecol Evol* 14:1065–1078.
- Davison, R. et al. (2019) Contributions of covariance. *Methods Ecol Evol* 10:1656–1672.
- Poelwijk, F.J. et al. (2016) Learning the pattern of epistasis. *PLoS Comput Biol* 12:e1004771.
