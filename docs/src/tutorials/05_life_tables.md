# Life Tables and Age-from-Stage Analysis

## Overview

Stage-structured MPMs classify individuals by developmental stage, not age. Yet many demographic quantities — survivorship $l(x)$, age-specific fecundity $m(x)$, life expectancy — are fundamentally age-based. The **age-from-stage** method (Cochran & Ellner 1992, Caswell 2001) extracts age-specific schedules from stage-classified models by tracking a cohort through the $\mathbf{U}$ matrix. This vignette demonstrates life table construction, survivorship conversions, and comparative survivorship analysis.

## Setup

```@example mpm
using MatrixProjectionModels
using Plots
using LinearAlgebra
```

## Example Species

We use the loggerhead sea turtle model — a classic long-lived species with delayed maturity.

```@example mpm
U_turtle = [0.0     0.0     0.0     0.0     0.0
            0.6747  0.7370  0.0     0.0     0.0
            0.0     0.0486  0.6610  0.0     0.0
            0.0     0.0     0.0147  0.6907  0.0
            0.0     0.0     0.0     0.0518  0.8091]

F_turtle = [0.0  0.0  0.0  0.0  127.0
            0.0  0.0  0.0  0.0  0.0
            0.0  0.0  0.0  0.0  0.0
            0.0  0.0  0.0  0.0  0.0
            0.0  0.0  0.0  0.0  0.0]

R_turtle = F_turtle  # No clonal reproduction
```

## Age-Specific Survivorship

The survivorship schedule $l(x)$ gives the probability of surviving from birth (stage 1) to age $x$. It is computed by tracking a cohort vector through repeated multiplication by $\mathbf{U}$:

$$l(x) = \mathbf{e}_1^\top \mathbf{U}^x \mathbf{e}_1$$

where the sum across stages gives the total surviving fraction at age $x$.

```@example mpm
lx = mpm_to_lx(U_turtle; start=1, xmax=80)

plot(0:length(lx)-1, lx,
    xlabel="Age (years)", ylabel="Survivorship l(x)",
    title="Survivorship — Loggerhead sea turtle",
    label=false, linewidth=2, color=:teal)
```

## Age-Specific Survival Probability

The per-age survival probability $p(x) = l(x+1) / l(x)$:

```@example mpm
px = mpm_to_px(U_turtle; start=1, xmax=80)

plot(0:length(px)-1, px,
    xlabel="Age (years)", ylabel="Survival probability p(x)",
    title="Age-specific survival — Sea turtle",
    label=false, linewidth=2, color=:blue)
```

## Age-Specific Hazard Rate

The hazard rate $h(x) = -\log(p(x))$ is the instantaneous mortality rate:

```@example mpm
hx = mpm_to_hx(U_turtle; start=1, xmax=80)

plot(0:length(hx)-1, hx,
    xlabel="Age (years)", ylabel="Hazard rate h(x)",
    title="Mortality hazard — Sea turtle",
    label=false, linewidth=2, color=:red)
```

## Age-Specific Reproduction

The reproductive schedule $m(x)$ gives the expected number of offspring per individual at age $x$:

```@example mpm
mx = mpm_to_mx(U_turtle, R_turtle; start=1, xmax=80)

plot(0:length(mx)-1, mx,
    xlabel="Age (years)", ylabel="Fecundity m(x)",
    title="Age-specific reproduction — Sea turtle",
    label=false, linewidth=2, color=:orange)
```

The delayed onset of reproduction reflects the long juvenile period — sea turtles take decades to reach maturity.

## Full Life Table

The `mpm_to_table` function returns all schedules in a single table:

```@example mpm
table = mpm_to_table(U_turtle, R_turtle; start=1, xmax=50)
for field in keys(table)
    println("$field: length = $(length(table[field]))")
end
```

```@example mpm
p1 = plot(table.x, table.lx, xlabel="Age", ylabel="l(x)",
    title="Survivorship", label=false, linewidth=2)
p2 = plot(table.x, table.px, xlabel="Age", ylabel="p(x)",
    title="Survival probability", label=false, linewidth=2, color=:blue)
p3 = plot(table.x, table.hx, xlabel="Age", ylabel="h(x)",
    title="Hazard rate", label=false, linewidth=2, color=:red)
p4 = plot(table.x, table.mx, xlabel="Age", ylabel="m(x)",
    title="Reproduction", label=false, linewidth=2, color=:orange)
plot(p1, p2, p3, p4, layout=(2,2), size=(800, 600))
```

## Life Table Conversions

The package provides bidirectional conversions between survivorship, survival probability, and hazard rate:

```@example mpm
# Start from survivorship
lx_orig = mpm_to_lx(U_turtle; xmax=50)

# Convert to survival probability and back
px_conv = lx_to_px(lx_orig)
lx_recov = px_to_lx(px_conv)
println("l(x) → p(x) → l(x) round-trip error: ",
    maximum(abs.(lx_recov[1:length(lx_orig)] - lx_orig)))

# Convert to hazard and back
hx_conv = lx_to_hx(lx_orig)
lx_from_h = hx_to_lx(hx_conv)
println("l(x) → h(x) → l(x) round-trip error: ",
    maximum(abs.(lx_from_h[1:length(lx_orig)] - lx_orig)))
```

## Comparative Survivorship

Survivorship curves reveal contrasting life history strategies. Deevey (1947) classified them into three types:

- **Type I** — low early mortality, high late-life mortality (convex curve; mammals)
- **Type II** — constant mortality throughout life (linear on log scale; some birds, reptiles)
- **Type III** — very high early mortality, low adult mortality (concave; many plants, fish, turtles)

```@example mpm
# Type I-like: large mammal (wolf-like from COMADRE)
U_mammal = [0.0   0.0   0.0   0.0
            0.50  0.0   0.0   0.0
            0.0   0.85  0.90  0.0
            0.0   0.0   0.08  0.80]

# Type II-like: constant-rate bird
U_bird = [0.0   0.0   0.0
          0.60  0.0   0.0
          0.0   0.70  0.70]

# Type III-like: sea turtle (from above)
U_type3 = U_turtle
```

```@example mpm
lx_mammal = mpm_to_lx(U_mammal; start=1, xmax=30)
lx_bird = mpm_to_lx(U_bird; start=1, xmax=30)
lx_turtle = mpm_to_lx(U_type3; start=1, xmax=80)

p = plot(xlabel="Relative age", ylabel="Survivorship l(x)",
    title="Deevey survivorship curves", legend=:topright, yscale=:log10,
    ylims=(1e-4, 1.5))

# Normalize ages to [0, 1] for comparison
ages_m = range(0, 1; length=length(lx_mammal))
ages_b = range(0, 1; length=length(lx_bird))
ages_t = range(0, 1; length=length(lx_turtle))

plot!(p, ages_m, lx_mammal, label="Type I (mammal)", linewidth=2)
plot!(p, ages_b, lx_bird, label="Type II (bird)", linewidth=2)
plot!(p, ages_t, lx_turtle, label="Type III (turtle)", linewidth=2)
p
```

## Starting from Different Stages

The `start` parameter specifies the initial stage for the age-from-stage calculation. This matters for species where individuals may enter the observed population at different stages:

```@example mpm
# Sea turtle: survivorship starting from different stages
p = plot(xlabel="Age", ylabel="l(x)",
    title="Survivorship by starting stage", legend=:topright)

stage_labels = ["Egg/hatchling", "Small juv", "Large juv", "Subadult", "Adult"]
for (i, label) in enumerate(stage_labels)
    lx_i = mpm_to_lx(U_turtle; start=i, xmax=40)
    plot!(p, 0:length(lx_i)-1, lx_i, label=label, linewidth=2)
end
p
```

Later stages have higher initial survival because they have already passed through the dangerous juvenile period.

## Summary

In this vignette we:

1. Extracted age-specific survivorship, survival probability, hazard rate, and fecundity from stage-structured MPMs
2. Built complete life tables via `mpm_to_table`
3. Demonstrated bidirectional conversions between $l(x)$, $p(x)$, and $h(x)$
4. Compared Deevey survivorship curve types across taxa
5. Explored the effect of starting stage on age-from-stage calculations

The next vignette covers life history traits derived from these schedules.
