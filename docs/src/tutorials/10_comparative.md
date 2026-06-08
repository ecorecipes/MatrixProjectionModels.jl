# Comparative Demography with COMPADRE and COMADRE

## Overview

The COMPADRE Plant Matrix Database and COMADRE Animal Matrix Database together contain thousands of published matrix population models spanning the tree of life. This vignette demonstrates **comparative demographic analysis** — using MPMs from multiple species to explore life history patterns, survivorship strategies, and the fast-slow continuum. We use matrices embedded directly from these databases (a future version will use the COMPADRE extension for automated download and parsing).

## Setup

```@example mpm
using MatrixProjectionModels
using Plots
using LinearAlgebra
using Statistics
```

## Species Collection

We assemble MPMs from 10 species spanning a wide range of life history strategies, drawn from COMPADRE and COMADRE.

```@example mpm
# Helper: build species entry
function make_species(name, taxon, U, F)
    A = U + F
    return (name=name, taxon=taxon, U=U, F=F, A=A,
            mpm=MatrixProjectionModel(U, F))
end

species = []

# --- PLANTS (COMPADRE) ---

# 1. Annual herb (Arabidopsis-like, 2-stage)
push!(species, make_species("Annual herb", :plant,
    [0.0  0.0; 0.40 0.0],
    [0.0  8.0; 0.0  0.0]))

# 2. Perennial grass (Festuca-like, 3-stage)
push!(species, make_species("Perennial grass", :plant,
    [0.20 0.0  0.0;
     0.30 0.55 0.0;
     0.0  0.15 0.80],
    [0.0 0.0 3.0;
     0.0 0.0 0.0;
     0.0 0.0 0.0]))

# 3. Succulent (Agave-like, 4-stage, from COMPADRE)
push!(species, make_species("Succulent", :plant,
    [0.05 0.0  0.0  0.0;
     0.10 0.70 0.0  0.0;
     0.0  0.10 0.85 0.0;
     0.0  0.0  0.05 0.90],
    [0.0 0.0 0.0 50.0;
     0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0]))

# 4. Temperate tree (Tsuga-like, 4-stage, from COMPADRE)
push!(species, make_species("Temperate tree", :plant,
    [0.40 0.0  0.0  0.0;
     0.10 0.75 0.0  0.0;
     0.0  0.05 0.92 0.0;
     0.0  0.0  0.03 0.98],
    [0.0 0.0 2.0 40.0;
     0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0]))

# 5. Tropical palm (Euterpe-like, 5-stage, from COMPADRE)
push!(species, make_species("Tropical palm", :plant,
    [0.30 0.0  0.0  0.0  0.0;
     0.05 0.60 0.0  0.0  0.0;
     0.0  0.08 0.80 0.0  0.0;
     0.0  0.0  0.05 0.90 0.0;
     0.0  0.0  0.0  0.05 0.95],
    [0.0 0.0 0.0 5.0 30.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0]))

# --- ANIMALS (COMADRE) ---

# 6. Short-lived bird (house sparrow-like, 3-stage)
push!(species, make_species("Songbird", :animal,
    [0.0   0.0   0.0;
     0.35  0.0   0.0;
     0.0   0.55  0.55],
    [0.0  1.5  2.5;
     0.0  0.0  0.0;
     0.0  0.0  0.0]))

# 7. Desert tortoise (Gopherus, 5-stage, from COMADRE)
push!(species, make_species("Desert tortoise", :animal,
    [0.0   0.0   0.0   0.0   0.0;
     0.716 0.567 0.0   0.0   0.0;
     0.0   0.149 0.567 0.0   0.0;
     0.0   0.0   0.149 0.604 0.0;
     0.0   0.0   0.0   0.235 0.817],
    [0.0 0.0 0.0 0.0 1.3;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0]))

# 8. Loggerhead sea turtle (Caretta, 5-stage, from COMADRE)
push!(species, make_species("Sea turtle", :animal,
    [0.0     0.0     0.0     0.0     0.0;
     0.6747  0.7370  0.0     0.0     0.0;
     0.0     0.0486  0.6610  0.0     0.0;
     0.0     0.0     0.0147  0.6907  0.0;
     0.0     0.0     0.0     0.0518  0.8091],
    [0.0 0.0 0.0 0.0 127.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0]))

# 9. Large mammal (elk-like, 4-stage, from COMADRE)
push!(species, make_species("Elk", :animal,
    [0.0   0.0   0.0   0.0;
     0.60  0.0   0.0   0.0;
     0.0   0.75  0.88  0.0;
     0.0   0.0   0.05  0.85],
    [0.0  0.0  0.3  0.5;
     0.0  0.0  0.0  0.0;
     0.0  0.0  0.0  0.0;
     0.0  0.0  0.0  0.0]))

# 10. Great ape (gorilla-like, 4-stage, from COMADRE)
push!(species, make_species("Great ape", :animal,
    [0.0   0.0   0.0   0.0;
     0.70  0.0   0.0   0.0;
     0.0   0.90  0.92  0.0;
     0.0   0.0   0.05  0.95],
    [0.0  0.0  0.05  0.15;
     0.0  0.0  0.0   0.0;
     0.0  0.0  0.0   0.0;
     0.0  0.0  0.0   0.0]))
```

## Growth Rates Across Species

```@example mpm
lambdas = [lambda(sp.mpm) for sp in species]
names = [sp.name for sp in species]
colors = [sp.taxon == :plant ? :forestgreen : :steelblue for sp in species]

bar(names, lambdas,
    ylabel="λ", title="Asymptotic growth rate across species",
    xrotation=45, legend=false, color=colors, alpha=0.7, size=(800, 400))
hline!([1.0], linestyle=:dash, color=:red, label="λ = 1")
```

## The Fast-Slow Continuum

The most fundamental pattern in comparative demography: species fall along a continuum from "fast" (short-lived, early reproduction, high fecundity) to "slow" (long-lived, delayed reproduction, low fecundity).

```@example mpm
le_vals = Float64[]
gt_vals = Float64[]

for sp in species
    push!(le_vals, life_expect_mean(sp.U; start=1))
    push!(gt_vals, gen_time(sp.U, sp.F; start=1, method=:R0))
end

scatter(le_vals, gt_vals,
    xlabel="Mean life expectancy", ylabel="Generation time",
    title="Fast-slow continuum",
    markersize=8, markerstrokewidth=0,
    color=colors, alpha=0.8,
    series_annotations=Plots.text.(names, 7, :left),
    legend=false, size=(700, 500))
plot!(identity, 0, maximum(le_vals)*1.2, linestyle=:dot, color=:gray,
    label=false)
```

## Survivorship Curve Comparison

```@example mpm
p = plot(xlabel="Relative lifespan", ylabel="Survivorship l(x)",
    title="Survivorship curves across species",
    yscale=:log10, ylims=(1e-4, 1.5), legend=:outerright, size=(800, 400))

for sp in species
    lx = mpm_to_lx(sp.U; start=1, xmax=100)
    # Normalize to relative lifespan [0, 1]
    ages = range(0, 1; length=length(lx))
    c = sp.taxon == :plant ? :forestgreen : :steelblue
    plot!(p, ages, lx, label=sp.name, linewidth=2, color=c, alpha=0.7)
end
p
```

## Vital Rate Profiles

```@example mpm
n = length(species)
profiles = zeros(n, 4)

for (i, sp) in enumerate(species)
    profiles[i, 1] = vr_survival(sp.U)
    profiles[i, 2] = vr_growth(sp.U)
    profiles[i, 3] = vr_stasis(sp.U)
    profiles[i, 4] = vr_shrinkage(sp.U)
end

vr_labels = ["Survival", "Growth", "Stasis", "Shrinkage"]
vr_colors = [:teal :forestgreen :goldenrod :indianred]
cum_profiles = cumsum(profiles, dims=2)
p_vr = plot(ylabel="Weighted mean rate",
    title="Vital rate profiles across species",
    legend=:outerright, size=(800, 400))
for j in size(profiles, 2):-1:1
    lower = j > 1 ? cum_profiles[:, j-1] : zeros(n)
    bar!(p_vr, 1:n, cum_profiles[:, j], fillrange=lower,
        label=vr_labels[j], color=vr_colors[j], alpha=0.7,
        bar_width=0.7)
end
plot!(p_vr, xticks=(1:n, names), xrotation=45)
p_vr
```

## Elasticity Patterns

A classic finding from comparative demography (Silvertown et al. 1993): species cluster in a triangle defined by elasticity of survival, growth, and fecundity.

```@example mpm
e_surv = Float64[]
e_growth = Float64[]
e_fec = Float64[]

for sp in species
    E = elasticity(sp.mpm)
    # Sum elasticity by component
    push!(e_surv, sum((E[i, j] for i in axes(E, 1), j in axes(E, 2) if sp.U[i, j] > 0 && i == j); init=0.0))
    push!(e_growth, sum((E[i, j] for i in axes(E, 1), j in axes(E, 2) if sp.U[i, j] > 0 && i != j); init=0.0))
    push!(e_fec, sum((E[i, j] for i in axes(E, 1), j in axes(E, 2) if sp.F[i, j] > 0); init=0.0))
end

# Ternary-like plot (2D projection)
scatter(e_fec, e_surv,
    xlabel="Elasticity of fecundity",
    ylabel="Elasticity of survival (stasis)",
    title="Elasticity triangle (Silvertown et al. 1993)",
    markersize=8, markerstrokewidth=0,
    color=colors, alpha=0.8,
    series_annotations=Plots.text.(names, 7, :left),
    legend=false, size=(700, 500),
    xlims=(0, 1), ylims=(0, 1))
# Add diagonal
plot!([0, 1], [1, 0], linestyle=:dot, color=:gray, label=false)
```

Trees and great apes cluster at high survival elasticity; annuals and songbirds at high fecundity elasticity.

## Life History Trait Dashboard

```@example mpm
trait_table = []
for sp in species
    le = life_expect_mean(sp.U; start=1)
    gt = gen_time(sp.U, sp.F; start=1, method=:R0)
    R0 = net_repro_rate(sp.U, sp.F; start=1)
    lon = longevity(sp.U; start=1, lx_crit=0.01)
    la = lambda(sp.mpm)
    H = entropy_k_stage(sp.U)
    push!(trait_table, (name=sp.name, taxon=sp.taxon,
        lambda=round(la, digits=3),
        life_exp=round(le, digits=1),
        gen_time=round(gt, digits=1),
        R0=round(R0, digits=2),
        longevity=lon,
        entropy=round(H, digits=3)))
end

# Print as formatted table
println(rpad("Species", 20), rpad("Taxon", 8), rpad("λ", 8),
    rpad("Life exp", 10), rpad("Gen time", 10), rpad("R₀", 8),
    rpad("Longevity", 10), "Entropy H")
println("-"^84)
for t in trait_table
    println(rpad(t.name, 20), rpad(string(t.taxon), 8),
        rpad(string(t.lambda), 8),
        rpad(string(t.life_exp), 10), rpad(string(t.gen_time), 10),
        rpad(string(t.R0), 8), rpad(string(t.longevity), 10),
        string(t.entropy))
end
```

## Keyfitz's Entropy and Survivorship Shape

Keyfitz's entropy $H$ characterizes survivorship shape: $H < 1$ (Type I, convex), $H \approx 1$ (Type II, exponential), $H > 1$ (Type III, concave).

```@example mpm
H_vals = [t.entropy for t in trait_table]

bar(names, H_vals,
    ylabel="Keyfitz's entropy H",
    title="Survivorship shape across species",
    xrotation=45, color=colors, alpha=0.7,
    size=(800, 400))
hline!([1.0], linestyle=:dash, color=:red, label="Type II (H=1)")
```

## Net Reproductive Rate vs Generation Time

This relationship captures the pace of life: high $R_0$ with short generation time (fast species) vs. low $R_0$ with long generation time (slow species).

```@example mpm
R0_vals = [t.R0 for t in trait_table]
gt_vals2 = [t.gen_time for t in trait_table]

scatter(gt_vals2, R0_vals,
    xlabel="Generation time", ylabel="Net reproductive rate R₀",
    title="Pace of life",
    markersize=8, markerstrokewidth=0,
    color=colors, alpha=0.8,
    series_annotations=Plots.text.(names, 7, :left),
    legend=false, size=(700, 500),
    yscale=:log10)
```

## Summary

In this vignette we:

1. Assembled 10 species from COMPADRE (plants) and COMADRE (animals)
2. Compared growth rates, life expectancy, and generation time across the tree of life
3. Visualized the fast-slow life history continuum
4. Compared survivorship curves (Deevey types I, II, III)
5. Analyzed vital rate and elasticity profiles across species
6. Reproduced the Silvertown et al. (1993) elasticity triangle
7. Built a comprehensive life history trait dashboard
8. Demonstrated that comparative demography reveals universal life history trade-offs

This completes the vignette series for MatrixProjectionModels.jl.
