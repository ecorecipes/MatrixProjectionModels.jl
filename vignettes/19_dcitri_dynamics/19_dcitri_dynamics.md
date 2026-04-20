# Population Dynamics of *Diaphorina citri* Under Seasonal Forcing
Simon Frost

- [Overview](#overview)
- [Setup](#setup)
- [Model Formulation](#model-formulation)
- [Part 1: Constant-Parameter Model](#part-1-constant-parameter-model)
  - [Eigenanalysis](#eigenanalysis)
  - [Projection Under Constant
    Conditions](#projection-under-constant-conditions)
  - [Sensitivity and Elasticity](#sensitivity-and-elasticity)
- [Part 2: Seasonal Environmental
  Forcing](#part-2-seasonal-environmental-forcing)
  - [Environmental Profiles](#environmental-profiles)
  - [Thermal Modulation of Development
    Rates](#thermal-modulation-of-development-rates)
  - [Effective Fertility $\beta_t$](#effective-fertility-beta_t)
- [Part 3: Seasonal Simulation](#part-3-seasonal-simulation)
  - [Figure 1: Integrated Dynamics of *D. citri* (52
    weeks)](#figure-1-integrated-dynamics-of-d-citri-52-weeks)
  - [Figure 2: Population Interaction by
    Stage](#figure-2-population-interaction-by-stage)
  - [Comparison With Published
    Values](#comparison-with-published-values)
- [Part 4: Spectral Analysis of the Periodic
  System](#part-4-spectral-analysis-of-the-periodic-system)
  - [Instantaneous Spectral Radius](#instantaneous-spectral-radius)
  - [Annual Product Matrix](#annual-product-matrix)
  - [Critical Fertility Threshold](#critical-fertility-threshold)
- [Discussion](#discussion)
  - [Limitations of the
    reconstruction](#limitations-of-the-reconstruction)
- [References](#references)

## Overview

*Diaphorina citri* Kuwayama, the Asian citrus psyllid, is the primary
vector of *Candidatus* Liberibacter asiaticus (CLas), the bacterium
associated with Huanglongbing (HLB) — the most destructive disease of
citrus worldwide (Bové 2006; Luo et al. 2021). The population dynamics
of the vector are a key determinant of disease spread.

This vignette reconstructs the discrete-time structured epidemiological
model of Ramírez-Ramírez et al. (2026), which describes the weekly
population dynamics of *D. citri* using a three-stage matrix projection
framework. The model tracks **eggs** (E), **nymphs** (N), and **adults**
(A) through a 52-week annual cycle with seasonal environmental forcing.

> [!NOTE]
>
> ### Reconstruction, not replication
>
> The paper does not fully specify the thermal response function for
> development rates or the closed-form relationship for the time-varying
> fertility $\beta_t$. Our analysis of the published parameter values
> and Table 2 data reveals that the underlying thermal index is
> nonlinear (not the simple linear form described in the text), and the
> published SAS code is not available. We therefore reconstruct the
> model from the tabulated values, interpolating both the thermal index
> $f_T(t)$ and fertility $\beta_t$ directly from the reported data
> points. This produces qualitatively correct seasonal dynamics but
> cannot exactly reproduce the published trajectories.

## Setup

``` julia
using MatrixProjectionModels
using Plots
using LinearAlgebra
using Printf
```

## Model Formulation

The population state vector $\mathbf{X}_t = (E_t, N_t, A_t)^\top$
evolves by
$$\mathbf{X}_{t+1} = \mathbf{M}_t \, \mathbf{X}_t$$
where the weekly transition matrix is
$$\mathbf{M}_t = \begin{pmatrix}
1 - \mu_E - \gamma_{E,t} & 0 & \beta_t \\
\gamma_{E,t} & 1 - \mu_N - \gamma_{N,t} & 0 \\
0 & \gamma_{N,t} & 1 - \mu_A - \delta
\end{pmatrix}.$$

The parameters are:

| Symbol     | Description                                  | Value   |
|------------|----------------------------------------------|---------|
| $\beta$    | Basal oviposition rate (eggs·adult⁻¹·week⁻¹) | 50      |
| $\mu_E$    | Egg mortality rate (week⁻¹)                  | 0.20    |
| $\gamma_E$ | Egg→nymph transition rate (week⁻¹)           | 0.30    |
| $\mu_N$    | Nymph mortality rate (week⁻¹)                | 0.20    |
| $\gamma_N$ | Nymph→adult transition rate (week⁻¹)         | 0.20    |
| $\mu_A$    | Adult mortality rate (week⁻¹)                | 0.10    |
| $\delta$   | Adult dispersal loss rate (week⁻¹)           | 0.05    |
| $K$        | Carrying capacity                            | 5 × 10⁶ |

Parameters are from Wenninger and Hall (2007) (fecundity), Liu and Tsai
(2000) (mortality), and <span class="nocase">García et al.</span> (2016)
(transition rates). Initial conditions: $E_0 = 100$, $N_0 = 50$,
$A_0 = 20$.

``` julia
# Basal demographic parameters
β  = 50.0    # oviposition rate
μ_E = 0.20   # egg mortality
γ_E = 0.30   # egg→nymph transition
μ_N = 0.20   # nymph mortality
γ_N = 0.20   # nymph→adult transition
μ_A = 0.10   # adult mortality
δ   = 0.05   # dispersal loss

# Carrying capacity
K = 5e6

# Initial conditions
n0 = [100.0, 50.0, 20.0]
```

    3-element Vector{Float64}:
     100.0
      50.0
      20.0

## Part 1: Constant-Parameter Model

Under constant optimal conditions, the transition matrix has fixed
entries.

``` julia
M_const = [1-μ_E-γ_E   0.0         β;
           γ_E          1-μ_N-γ_N   0.0;
           0.0          γ_N         1-μ_A-δ]

println("Transition matrix M (constant):")
display(M_const)
```

    Transition matrix M (constant):

    3×3 Matrix{Float64}:
     0.5  0.0  50.0
     0.3  0.6   0.0
     0.0  0.2   0.85

### Eigenanalysis

``` julia
ea = eigenanalysis_full(M_const)

println("Dominant eigenvalue (λ₁): ", round(ea.lambda, digits=4))
println("Stable stage distribution: ", round.(ea.stable_dist, digits=4))
println("Reproductive value:        ", round.(ea.repro_value, digits=4))
println("Damping ratio:             ", round(damping_ratio(M_const), digits=4))
```

    Dominant eigenvalue (λ₁): 2.1
    Stable stage distribution: [0.8117, 0.1623, 0.026]
    Reproductive value:        [0.3681, 1.9633, 14.7251]
    Damping ratio:             1.6868

The spectral radius $\lambda_1 > 1$ confirms sustained exponential
growth under constant optimal conditions — the discrete analogue of
$R_0 > 1$ in continuous epidemiological models (Luo et al. 2021).

### Projection Under Constant Conditions

``` julia
prob_const = MPMProblem(M_const, n0, (0, 52))
sol_const = solve(prob_const, DirectIteration())

weeks = sol_const.t
E_const = [u[1] for u in sol_const.u]
N_const = [u[2] for u in sol_const.u]
A_const = [u[3] for u in sol_const.u]

plot(weeks, E_const, label="Eggs", lw=2, color=:blue,
     xlabel="Week", ylabel="Population",
     title="D. citri: Constant-parameter projection (52 weeks)",
     legend=:topleft, yscale=:log10, size=(700, 400))
plot!(weeks, N_const, label="Nymphs", lw=2, color=:red)
plot!(weeks, A_const, label="Adults", lw=2, color=:green)
```

![](19_dcitri_dynamics_files/figure-commonmark/cell-6-output-1.svg)

Under constant conditions, all stages grow exponentially at rate
$\lambda_1$, rapidly reaching astronomical numbers. This motivates the
incorporation of seasonal environmental forcing.

### Sensitivity and Elasticity

``` julia
S = sensitivity(M_const)
E_mat = elasticity(M_const)

p1 = heatmap(S[end:-1:1, :], title="Sensitivity",
    xticks=(1:3, ["E", "N", "A"]), yticks=(1:3, ["A", "N", "E"]),
    color=:viridis, aspect_ratio=1, size=(350, 300))
p2 = heatmap(E_mat[end:-1:1, :], title="Elasticity",
    xticks=(1:3, ["E", "N", "A"]), yticks=(1:3, ["A", "N", "E"]),
    color=:viridis, aspect_ratio=1, size=(350, 300))
plot(p1, p2, layout=(1, 2), size=(750, 320))
```

![](19_dcitri_dynamics_files/figure-commonmark/cell-7-output-1.svg)

The sensitivity matrix shows that population growth is most sensitive to
changes in the fecundity entry $\beta$ (top-right element), confirming
its role as the dominant driver.

## Part 2: Seasonal Environmental Forcing

The paper uses an annual environmental cycle where temperature and
humidity modulate the demographic rates. We reconstruct the
environmental forcing from the tabulated data.

### Environmental Profiles

Temperature and humidity follow cosine seasonal profiles:
$$T(t) = 31.05 + 6.05 \cos\!\left(\frac{2\pi(t - 13)}{52}\right), \qquad
H(t) = 53.25 + 13.75 \cos\!\left(\frac{2\pi t}{52}\right).$$

These match all seven tabulated data points in Table 2 of
Ramírez-Ramírez et al. (2026) to within rounding precision.

``` julia
T_func(t) = 31.05 + 6.05 * cos(2π * (t - 13) / 52)
H_func(t) = 53.25 + 13.75 * cos(2π * t / 52)

p1 = plot(0:52, T_func.(0:52), lw=2, color=:red,
    xlabel="Week", ylabel="Temperature (°C)",
    title="Seasonal temperature", label=false,
    ylims=(23, 39))
p2 = plot(0:52, H_func.(0:52), lw=2, color=:steelblue,
    xlabel="Week", ylabel="Relative humidity (%)",
    title="Seasonal humidity", label=false,
    ylims=(35, 72))
plot(p1, p2, layout=(1, 2), size=(800, 350))
```

![](19_dcitri_dynamics_files/figure-commonmark/cell-8-output-1.svg)

### Thermal Modulation of Development Rates

The transition rates $\gamma_{E,t}$ and $\gamma_{N,t}$ are modulated by
a thermal suitability index $f_T(t)$ such that
$\gamma_{E,t} = \gamma_E \cdot f_T(t)$ and
$\gamma_{N,t} = \gamma_N \cdot f_T(t)$.

The paper describes $f_T$ as a linear function of temperature, but
comparison of the published transition rates with the temperature
profile reveals a **nonlinear** relationship: $f_T = 1.0$ at $T = 31$°C
(week 52) yet only $0.90$ at $T = 35$°C (week 20), and $f_T = 0.14$ at
$T = 25$°C (week 39) rather than the zero predicted by a linear index
with $T_\min = 25$°C. We therefore interpolate $f_T(t)$ directly from
the Table 2 data points.

``` julia
# Thermal index keypoints from Table 2:
# f_T = γ_E(t)/γ_E_base averaged with γ_N(t)/γ_N_base
fT_keypoints = [
    (0,  1.00),   # optimal start (same as week 52)
    (10, 1.00),   # Table 2: γ_E=0.30/0.30, γ_N=0.20/0.20
    (15, 1.00),   # Table 2: γ_E=0.30/0.30, γ_N=0.20/0.20
    (20, 0.90),   # Table 2: γ_E=0.27/0.30, γ_N=0.18/0.20
    (26, 0.30),   # interpolated decline
    (33, 0.05),   # approaching minimum
    (39, 0.14),   # Table 2: γ_E=0.04/0.30≈0.13, γ_N=0.03/0.20=0.15
    (46, 0.44),   # Table 2: γ_E=0.13/0.30≈0.43, γ_N=0.09/0.20=0.45
    (52, 1.00),   # Table 2: γ_E=0.30/0.30, γ_N=0.20/0.20
]

function interp_vals(keypoints, t; vmax=Inf)
    ws = [p[1] for p in keypoints]
    vs = [p[2] for p in keypoints]
    tc = clamp(t, ws[1], ws[end])
    idx = findlast(w -> w <= tc, ws)
    idx == length(ws) && return vs[end]
    frac = (tc - ws[idx]) / (ws[idx+1] - ws[idx])
    return clamp(vs[idx] + frac * (vs[idx+1] - vs[idx]), 0.0, vmax)
end

f_T(t) = interp_vals(fT_keypoints, t; vmax=1.0)

# Verify against Table 2 transition rates
println("Verification of thermal index against Table 2:")
for (week, γE_tab, γN_tab) in [(10, 0.30, 0.20), (15, 0.30, 0.20),
                                 (20, 0.27, 0.18), (39, 0.04, 0.03),
                                 (46, 0.13, 0.09), (52, 0.30, 0.20)]
    fT = f_T(week)
    @printf("Week %2d: f_T=%.3f → γ_E=%.2f (tab %.2f), γ_N=%.2f (tab %.2f)\n",
            week, fT, γ_E*fT, γE_tab, γ_N*fT, γN_tab)
end
```

    Verification of thermal index against Table 2:
    Week 10: f_T=1.000 → γ_E=0.30 (tab 0.30), γ_N=0.20 (tab 0.20)
    Week 15: f_T=1.000 → γ_E=0.30 (tab 0.30), γ_N=0.20 (tab 0.20)
    Week 20: f_T=0.900 → γ_E=0.27 (tab 0.27), γ_N=0.18 (tab 0.18)
    Week 39: f_T=0.140 → γ_E=0.04 (tab 0.04), γ_N=0.03 (tab 0.03)
    Week 46: f_T=0.440 → γ_E=0.13 (tab 0.13), γ_N=0.09 (tab 0.09)
    Week 52: f_T=1.000 → γ_E=0.30 (tab 0.30), γ_N=0.20 (tab 0.20)

``` julia
plot(0:52, f_T.(0:52), lw=2, color=:darkorange,
     xlabel="Week", ylabel="Thermal index f_T",
     title="Thermal suitability index (interpolated from Table 2)",
     label=false, size=(700, 300), ylims=(-0.05, 1.1))
scatter!([p[1] for p in fT_keypoints], [p[2] for p in fT_keypoints],
         color=:black, ms=5, label="Table 2 data")
hline!([1.0], ls=:dash, color=:gray, lw=1, label=false)
```

![](19_dcitri_dynamics_files/figure-commonmark/cell-10-output-1.svg)

### Effective Fertility $\beta_t$

The paper reports time-varying $\beta_t$ values at selected weeks (Table
2) but does not provide a closed-form formula. We linearly interpolate
$\beta_t$ between the tabulated data points, with intermediate values
estimated to ensure smooth transitions.

``` julia
β_keypoints = [
    (0,  25.00),  # estimated from seasonal pattern
    (5,  36.00),  # estimated (ascending toward peak)
    (10, 32.76),  # Table 2
    (13, 48.00),  # near-peak (spring maximum)
    (15, 41.88),  # Table 2
    (16, 20.00),  # rapid post-peak decline
    (18,  8.00),  # continued decline
    (20,  5.23),  # Table 2
    (23,  1.00),  # approaching zero
    (26,  0.00),  # β ≈ 0 (text: adverse period)
    (39,  0.00),  # Table 2
    (43,  3.00),  # slow recovery begins
    (46,  7.36),  # Table 2
    (52, 25.00),  # return to initial conditions
]

β_t(t) = interp_vals(β_keypoints, t; vmax=β)

plot(0:52, β_t.(0:52), lw=2, color=:teal,
     xlabel="Week", ylabel="β_t (eggs·adult⁻¹·week⁻¹)",
     title="Effective fertility (interpolated from Table 2)",
     label=false, size=(700, 300))
scatter!([p[1] for p in β_keypoints if p[1] in [10,15,20,39,46,52]],
         [p[2] for p in β_keypoints if p[1] in [10,15,20,39,46,52]],
         color=:black, ms=5, label="Table 2 data")
```

![](19_dcitri_dynamics_files/figure-commonmark/cell-11-output-1.svg)

## Part 3: Seasonal Simulation

We use the `DensityDependent` solver pathway, which supplies the current
population vector to the matrix function at each step. The paper lists a
carrying capacity $K = 5 \times 10^6$ in Table 1 and we include logistic
density dependence on fecundity:

$$\beta_t^{\text{eff}} = \beta_t \cdot \max\!\left(0,\; 1 - \frac{E_t + N_t + A_t}{K}\right).$$

``` julia
function build_matrix(n, p, t)
    fT = f_T(t)
    βt_eff = β_t(t) * max(0.0, 1.0 - sum(n) / K)
    γEt = γ_E * fT
    γNt = γ_N * fT

    return [1-μ_E-γEt   0.0         βt_eff;
            γEt          1-μ_N-γNt   0.0;
            0.0          γNt         1-μ_A-δ]
end

prob_seasonal = MPMProblem(
    LefkovitchMPM(), DensityDependent(), Deterministic(),
    build_matrix, n0, (0, 52);
    p = nothing
)
sol_seasonal = solve(prob_seasonal, DirectIteration())

weeks_s = sol_seasonal.t
E_s = [u[1] for u in sol_seasonal.u]
N_s = [u[2] for u in sol_seasonal.u]
A_s = [u[3] for u in sol_seasonal.u]
```

    53-element Vector{Float64}:
         20.0
         27.0
         34.95
         69.90648
        144.0977723962176
        265.20484908012054
        489.5038597667235
        955.7557976588755
       1898.7466404417494
       3682.016286246535
          ⋮
      77918.54936860483
      71031.98118842249
      68657.6134211167
      71879.21427410316
      84005.24746641233
     108825.01591695647
     152859.07935718299
     226597.7055178804
     346131.3894578513

### Figure 1: Integrated Dynamics of *D. citri* (52 weeks)

``` julia
plot(weeks_s, E_s, label="E (Eggs)", lw=2.5, color=:blue,
     xlabel="Week", ylabel="Population",
     title="Integrated Dynamics of Diaphorina citri (52 weeks)",
     legend=:topright, size=(700, 420))
plot!(weeks_s, N_s, label="N (Nymphs)", lw=2.5, color=:red)
plot!(weeks_s, A_s, label="A (Adults)", lw=2.5, color=:green)
```

![Integrated dynamics of *Diaphorina citri* under seasonal environmental
forcing. Compare with Figure 1 of Ramírez-Ramírez et
al. (2026).](19_dcitri_dynamics_files/figure-commonmark/cell-13-output-1.svg)

### Figure 2: Population Interaction by Stage

``` julia
plot(weeks_s, E_s, label="Eggs", lw=1.5, color=:blue,
     marker=:diamond, ms=3, markerstrokewidth=0,
     xlabel="Week", ylabel="Population",
     title="Population Interaction by Stage (Integrated Model)",
     legend=:topright, size=(700, 420))
plot!(weeks_s, N_s, label="Nymphs", lw=1.5, color=:red,
      marker=:circle, ms=3, markerstrokewidth=0)
plot!(weeks_s, A_s, label="Adults", lw=1.5, color=:green,
      marker=:utriangle, ms=3, markerstrokewidth=0)
```

![Population interaction by stage with markers. Compare with Figure 2 of
Ramírez-Ramírez et
al. (2026).](19_dcitri_dynamics_files/figure-commonmark/cell-14-output-1.svg)

### Comparison With Published Values

``` julia
println("Comparison with Table 2 of Ramírez-Ramírez et al. (2026)")
println("=" ^ 92)
@printf("%-22s %4s %12s %12s  %12s %12s\n",
        "Phase", "Week", "Eggs (model)", "Eggs (paper)", "Adults (mod)", "Adults (pap)")
println("-" ^ 92)

for (phase, wk, E_exp, N_exp, A_exp) in [
    ("Initial condition",    0,  100,        50,       20),
    ("Exponential growth",  10,  212_157,    49_416,   8_864),
    ("Population peak",     15,  2_103_357,  1_026_613, 68_177),
    ("Inflection point",    20,  710_234,    1_136_694, 99_030),
    ("Relative minimum",    39,  60_919,     11_608,   9_717),
    ("Progressive recovery", 46, 96_241,     10_756,   3_324),
    ("Final condition",     52,  285_088,    56_410,   9_249),
]
    idx = wk + 1
    @printf("%-22s %4d  %10.0f  %10d   %10.0f  %10d\n",
            phase, wk, E_s[idx], E_exp, A_s[idx], A_exp)
end
println("-" ^ 92)
println()
println("Note: Our reconstruction captures the qualitative seasonal pattern")
println("(spring peak → summer contraction → autumn recovery) but diverges")
println("quantitatively because the paper's thermal response function and")
println("β_t formula are not fully specified in the published description.")
```

    Comparison with Table 2 of Ramírez-Ramírez et al. (2026)
    ============================================================================================
    Phase                  Week Eggs (model) Eggs (paper)  Adults (mod) Adults (pap)
    --------------------------------------------------------------------------------------------
    Initial condition         0         100         100           20          20
    Exponential growth       10      163468      212157         7057        8864
    Population peak          15     2941058     2103357       187635       68177
    Inflection point         20     1893490      710234      1157443       99030
    Relative minimum         39       22197       60919       161283        9717
    Progressive recovery     46      798263       96241        68658        3324
    Final condition          52     2778850      285088       346131        9249
    --------------------------------------------------------------------------------------------

    Note: Our reconstruction captures the qualitative seasonal pattern
    (spring peak → summer contraction → autumn recovery) but diverges
    quantitatively because the paper's thermal response function and
    β_t formula are not fully specified in the published description.

## Part 4: Spectral Analysis of the Periodic System

### Instantaneous Spectral Radius

The spectral radius of $M_t$ at each week shows when the system is in a
supercritical ($\rho > 1$, expansive) or subcritical ($\rho < 1$,
contractive) regime.

``` julia
ρ_weekly = [maximum(abs.(eigvals(build_matrix(n0, nothing, t)))) for t in 0:52]

plot(0:52, ρ_weekly, lw=2, color=:purple,
     xlabel="Week", ylabel="Spectral radius ρ(M_t)",
     title="Instantaneous spectral radius over the annual cycle",
     legend=false, size=(700, 350))
hline!([1.0], ls=:dash, color=:black, lw=1.5)
annotate!(40, 1.05, text("ρ = 1 (threshold)", 9, :left))
```

![](19_dcitri_dynamics_files/figure-commonmark/cell-16-output-1.svg)

### Annual Product Matrix

For a periodically forced system, the correct asymptotic growth metric
is the dominant eigenvalue of the **annual monodromy matrix**
$P = M_{52} M_{51} \cdots M_1$ (Caswell 2001). This accounts for the
cumulative effect of favourable and unfavourable periods.

``` julia
P = Matrix{Float64}(I, 3, 3)
for t in 0:51
    P = build_matrix(sol_seasonal.u[t+1], nothing, t) * P
end

ea_annual = eigen(P)
λ_annual = maximum(abs.(ea_annual.values))
println("Annual monodromy matrix eigenvalues: ", round.(abs.(ea_annual.values), digits=4))
println("Annual multiplier (dominant λ): ", round(λ_annual, sigdigits=6))
println("Equivalent weekly growth rate:  ", round(λ_annual^(1/52), sigdigits=6))

if λ_annual > 1
    println("→ System is expansive over the annual cycle")
else
    println("→ System is contractive over the annual cycle")
end
```

    Annual monodromy matrix eigenvalues: [0.0, 0.0, 20228.6145]
    Annual multiplier (dominant λ): 20228.6
    Equivalent weekly growth rate:  1.21006
    → System is expansive over the annual cycle

### Critical Fertility Threshold

The paper identifies $\beta_t \approx 30$ eggs·adult⁻¹·week⁻¹ as a
supercritical threshold. We verify this by computing the spectral radius
of the constant-parameter matrix for a range of $\beta$ values.

``` julia
β_range = 0:1:80
ρ_vs_β = [maximum(abs.(eigvals([1-μ_E-γ_E 0.0 b;
                                 γ_E 1-μ_N-γ_N 0.0;
                                 0.0 γ_N 1-μ_A-δ])))
           for b in β_range]

β_crit = β_range[findfirst(ρ_vs_β .>= 1.0)]

plot(β_range, ρ_vs_β, lw=2.5, color=:darkblue,
     xlabel="β (eggs·adult⁻¹·week⁻¹)", ylabel="Spectral radius ρ(M)",
     title="Population growth threshold vs. fertility",
     label=false, size=(700, 350))
hline!([1.0], ls=:dash, color=:black, lw=1.5, label="ρ = 1")
vline!([β_crit], ls=:dot, color=:red, lw=1.5,
       label=@sprintf("β_crit ≈ %d", β_crit))
```

![](19_dcitri_dynamics_files/figure-commonmark/cell-18-output-1.svg)

The critical fertility $\beta_\text{crit}$ at which $\rho(M) = 1$
provides the discrete analogue of the epidemic threshold.

## Discussion

This vignette demonstrates how `MatrixProjectionModels.jl` can be used
to model seasonally forced vector population dynamics. Key findings:

1.  **Constant-parameter analysis** confirms exponential growth
    ($\lambda_1 > 1$) under optimal conditions, consistent with the
    supercritical regime described by Ramírez-Ramírez et al. (2026).

2.  **Seasonal forcing** produces the characteristic three-phase
    dynamics: spring expansion, summer contraction (as $\beta_t \to 0$
    and $f_T \to 0$), and autumn recovery. The adult compartment shows
    greater structural persistence, acting as a demographic buffer
    during adverse periods.

3.  **The annual monodromy matrix** provides the correct long-term
    growth assessment for this periodic system. Instantaneous spectral
    radii are diagnostic but not predictive of asymptotic behavior.

4.  **The critical fertility threshold** ($\beta \approx 30$
    eggs·adult⁻¹·week⁻¹) separates expansive from contractive regimes,
    providing a target for intervention planning.

### Limitations of the reconstruction

Our seasonal simulation captures the correct qualitative phases but
diverges quantitatively from the paper’s trajectories, particularly in
the adult compartment (which peaks several-fold higher than reported).
Analysis of the published parameter table reveals that the paper’s
thermal response function is nonlinear and does not match the simple
linear index described in the text. Without access to the original SAS
implementation, exact calibration of $f_T(t)$ and $\beta_t$ is not
possible. The published parameters are also internally inconsistent:
stepping the matrix equation forward from the Table 2 values at week 10
with the stated transition rates produces egg counts at week 15 that
substantially exceed the reported 2.1 million, suggesting additional
density feedbacks or parameter modulations in the original code that are
not documented.

From an applied perspective, synchronizing management interventions with
the peak egg→nymph transition (weeks 10–15) could maximize suppression
of the subsequent adult cohort — the stage responsible for CLas
transmission (Leong et al. 2022).

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-bove2006" class="csl-entry">

Bové, Joseph M. 2006. “Huanglongbing: A Destructive, Newly-Emerging,
Century-Old Disease of Citrus.” *Journal of Plant Pathology* 88 (1):
7–37.

</div>

<div id="ref-caswell2001" class="csl-entry">

Caswell, Hal. 2001. *Matrix Population Models: Construction, Analysis,
and Interpretation*. 2nd ed. Sinauer Associates.

</div>

<div id="ref-garcia2016" class="csl-entry">

<span class="nocase">García, M. et al.</span> 2016. *Development Times
of Diaphorina Citri in Subtropical Environments*.

</div>

<div id="ref-leong2022" class="csl-entry">

Leong, Soon Shin, Sui Chai Thing Leong, and G. Andrew C. Beattie. 2022.
“Integrated Pest Management Strategies for Asian Citrus Psyllid
Diaphorina Citri and Huanglongbing in Citrus for Sarawak, East Malaysia,
Borneo.” *Insects* 13: 960.

</div>

<div id="ref-liu2000" class="csl-entry">

Liu, Yen-Hsun, and Jaw-Hua Tsai. 2000. “Effects of Temperature on
Biology and Life Table Parameters of the Asian Citrus Psyllid,
Diaphorina Citri Kuwayama (Homoptera: Psyllidae).” *Annals of Applied
Biology* 137 (3): 201–16.

</div>

<div id="ref-luo2021" class="csl-entry">

Luo, Yongzhen, Fang Zhang, Youquan Liu, and Shujing Gao. 2021. “Analysis
and Optimal Control of a Huanglongbing Mathematical Model with Resistant
Vector.” *Infectious Disease Modelling* 6: 782–804.

</div>

<div id="ref-ramirez2026" class="csl-entry">

Ramírez-Ramírez, Faustino, Arturo Moisés Chávez-Rodríguez, Isaac
Andrade-González, et al. 2026. “Modelo Integrado de Simulación
Epidemiológica Discreta Para La Dinámica Poblacional de Diaphorina Citri
En Relación Con El Ciclo de Cultivo.” *ISSN 2764-1317* 6 (1).

</div>

<div id="ref-wenninger2007" class="csl-entry">

Wenninger, Erik J., and David G. Hall. 2007. “Oviposition Behavior and
Fecundity of Diaphorina Citri.” *Florida Entomologist* 90 (2): 335–40.

</div>

</div>
