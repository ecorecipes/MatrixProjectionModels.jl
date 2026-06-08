# Stochastic Matrix Projection Models

## Overview

Real populations experience **environmental stochasticity** — year-to-year variation in vital rates due to weather, resource availability, or disturbance. Stochastic MPMs model this by drawing from a set of environment-specific matrices at each time step. The stochastic growth rate $\lambda_s \leq \lambda_{\text{det}}$ (Tuljapurkar's inequality), meaning environmental variation generally reduces long-term population growth. This vignette demonstrates both kernel-resampled and parameter-resampled stochastic MPMs.

## Setup

```@example mpm
using MatrixProjectionModels
using Plots
using LinearAlgebra
using Statistics
using Random
Random.seed!(42)
```

## The Model Species

We model a perennial plant (inspired by *Silene spaldingii*, Spalding's catchfly, from COMPADRE) with 4 stages: seedling, small, large, flowering. The species experiences alternating good and bad years driven by rainfall.

### Environment-Specific Matrices

```@example mpm
# Good year: high survival, high growth, high fecundity
U_good = [0.10  0.00  0.00  0.00
          0.25  0.50  0.00  0.00
          0.00  0.25  0.65  0.10
          0.00  0.00  0.20  0.80]
F_good = [0.0  0.0  0.0  60.0
          0.0  0.0  0.0  0.0
          0.0  0.0  0.0  0.0
          0.0  0.0  0.0  0.0]

# Bad year: low survival, reduced growth, low fecundity
U_bad = [0.05  0.00  0.00  0.00
         0.15  0.35  0.00  0.00
         0.00  0.10  0.45  0.15
         0.00  0.00  0.10  0.55]
F_bad = [0.0  0.0  0.0  15.0
         0.0  0.0  0.0  0.0
         0.0  0.0  0.0  0.0
         0.0  0.0  0.0  0.0]

A_good = U_good + F_good
A_bad = U_bad + F_bad

mpm_good = MatrixProjectionModel(U_good, F_good)
mpm_bad = MatrixProjectionModel(U_bad, F_bad)

println("Good year λ = ", round(lambda(mpm_good), digits=4))
println("Bad year λ  = ", round(lambda(mpm_bad), digits=4))
```

### Mean Matrix

The mean matrix $\bar{\mathbf{A}} = \frac{1}{2}(\mathbf{A}_{\text{good}} + \mathbf{A}_{\text{bad}})$ gives the deterministic growth rate if we ignored environmental variation:

```@example mpm
A_mean = 0.5 * A_good + 0.5 * A_bad
lambda_det = lambda(A_mean)
println("Deterministic λ (mean matrix) = ", round(lambda_det, digits=4))
```

## Kernel-Resampled Stochastic MPMs

In the **kernel-resampled** approach, we pre-compute matrices for each environment and randomly draw one at each time step.

```@example mpm
n0 = [100.0, 50.0, 25.0, 10.0]

prob_stoch = MPMProblem(
    StochasticKernelResampled(),
    [mpm_good, mpm_bad],
    n0, (0, 200))

sol_stoch = solve(prob_stoch, DirectIteration())
```

```@example mpm
total_pop = [sum(u) for u in sol_stoch.u]

plot(0:200, log10.(total_pop),
    xlabel="Time step (years)", ylabel="log₁₀ N(t)",
    title="Stochastic population trajectory",
    label="Stochastic", linewidth=1, alpha=0.7, color=:blue)

# Compare with deterministic trajectory
prob_det = MPMProblem(MatrixProjectionModel(A_mean), n0, (0, 200))
sol_det = solve(prob_det, DirectIteration())
total_det = [sum(u) for u in sol_det.u]
plot!(0:200, log10.(total_det),
    label="Deterministic (mean)", linewidth=2, color=:red, linestyle=:dash)
```

### Stochastic Growth Rate

The stochastic growth rate $\lambda_s$ is the geometric mean of per-step growth rates:

$$\log \lambda_s = \lim_{T \to \infty} \frac{1}{T} \sum_{t=1}^{T} \log \lambda(t)$$

We compute it from the simulation solution:

```@example mpm
lambda_s = stochastic_growth_rate(sol_stoch)
println("Stochastic growth rate (λₛ) = ", round(lambda_s, digits=4))
println("Deterministic growth rate (λ̄) = ", round(lambda_det, digits=4))
println("Tuljapurkar's inequality: λₛ ≤ λ̄ → ", lambda_s <= lambda_det + 1e-10)
```

### Multiple Realizations

```@example mpm
p = plot(xlabel="Time step", ylabel="log₁₀ N(t)",
    title="20 stochastic realizations", alpha=0.3)

for _ in 1:20
    sol_i = solve(prob_stoch, DirectIteration())
    pop_i = [sum(u) for u in sol_i.u]
    plot!(p, 0:200, log10.(pop_i), label=false, color=:blue, alpha=0.2)
end

# Overlay deterministic
plot!(p, 0:200, log10.(total_det),
    label="Deterministic", linewidth=2, color=:red, linestyle=:dash)
p
```

## Parameter-Resampled Stochastic MPMs

In the **parameter-resampled** approach, vital rate parameters are drawn from distributions at each time step. This allows continuous environmental variation rather than discrete states.

```@example mpm
# Environmental state function: env_state(rng, t) → parameter tuple
function env_resample(rng, t)
    # Rainfall modifier (log-normal: mean 1, CV 0.3)
    rain = exp(randn(rng) * 0.3 - 0.045)  # mean-corrected
    return (rain=rain,)
end

# Matrix function: matrix(params) → projection matrix
function matrix_fn(params)
    r = params.rain
    # Scale survival and fecundity by rainfall
    U_t = clamp.(U_good .* (0.5 + 0.5 * r), 0, 1)
    F_t = F_good .* r
    return U_t + F_t
end

prob_param = MPMProblem(
    StochasticParameterResampled(),
    matrix_fn,
    n0, (0, 200);
    env_state=env_resample)

sol_param = solve(prob_param, DirectIteration())
```

```@example mpm
total_param = [sum(u) for u in sol_param.u]

plot(0:200, log10.(total_param),
    xlabel="Time step", ylabel="log₁₀ N(t)",
    title="Parameter-resampled stochastic dynamics",
    label="Parameter resampled", linewidth=1, color=:purple)
plot!(0:200, log10.(total_det),
    label="Deterministic", linewidth=2, color=:red, linestyle=:dash)
```

## Stochastic Perturbation Analysis

We can compute the sensitivity and elasticity of the *stochastic* growth rate $\log \lambda_s$ to changes in mean matrix elements:

```@example mpm
stoch_sens = perturb_stochastic([A_good, A_bad]; type=:sensitivity, n_sim=50000)

stage_labels = ["Seedling", "Small", "Large", "Flowering"]
heatmap(stage_labels, stage_labels, stoch_sens,
    title="Stochastic sensitivity of log λₛ",
    xlabel="From stage", ylabel="To stage",
    color=:YlOrRd)
```

### Comparing Deterministic vs Stochastic Sensitivity

```@example mpm
S_det = sensitivity(MatrixProjectionModel(A_mean))

p1 = heatmap(stage_labels, stage_labels, S_det,
    title="Deterministic sensitivity", color=:YlOrRd)
p2 = heatmap(stage_labels, stage_labels, stoch_sens,
    title="Stochastic sensitivity", color=:YlOrRd)
plot(p1, p2, layout=(1,2), size=(800, 350))
```

Environmental stochasticity can shift sensitivity patterns — conservation priorities derived from deterministic models may not hold in variable environments.

## Effect of Variance on Growth

Tuljapurkar's approximation shows that $\log \lambda_s \approx \log \bar{\lambda} - \tau^2/(2\bar{\lambda}^2)$, where $\tau^2$ is the variance in $\lambda$ across environments. More variable environments reduce population growth more:

```@example mpm
# Vary the severity of bad years
bad_scalers = 0.5:0.1:1.0  # 1.0 = same as good year (no variation)
lambda_s_values = Float64[]

for scale in bad_scalers
    A_bad_scaled = (1 - scale) * A_bad + scale * A_good
    mpm_bad_s = MatrixProjectionModel(A_bad_scaled)
    prob_s = MPMProblem(
        StochasticKernelResampled(),
        [mpm_good, mpm_bad_s],
        n0, (0, 500))
    sol_s = solve(prob_s, DirectIteration())
    ls = stochastic_growth_rate(sol_s)
    push!(lambda_s_values, ls)
end

plot(1.0 .- bad_scalers, lambda_s_values,
    xlabel="Environmental variation (1 - scale)",
    ylabel="λₛ",
    title="Stochastic growth rate vs environmental variation",
    label=false, linewidth=2, color=:blue)
hline!([lambda(mpm_good)], label="λ (good year)", linestyle=:dash, color=:green)
```

## Summary

In this vignette we:

1. Modeled environmental stochasticity using good-year/bad-year matrices
2. Demonstrated kernel-resampled stochastic MPMs
3. Computed the stochastic growth rate $\lambda_s$ and verified Tuljapurkar's inequality
4. Visualized multiple stochastic realizations
5. Implemented parameter-resampled stochasticity with continuous environmental variation
6. Performed stochastic perturbation analysis
7. Showed how environmental variance reduces long-term population growth

The next vignette covers population simulation and density dependence.
