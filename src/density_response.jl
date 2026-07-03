"""
Density response functions for population projection models.

Provides standard functional forms for density-dependent vital rate modification:
- Ricker (overcompensatory)
- Beverton-Holt (compensatory)
- Usher/logistic (sigmoidal)
- Logistic (linear density dependence)
- Theta-logistic (generalized)

Each function takes population density N and returns a multiplicative modifier
for the corresponding vital rate.

Reference: lefko3 package (Shefferson & Ehrlen); Caswell (2001) Ch. 16.
"""

# --- Abstract type ---

"""
    AbstractDensityResponse

Abstract supertype for density-dependent response functions.
All subtypes must implement `(f::AbstractDensityResponse)(N::Real)` returning
a multiplicative modifier in [0, ∞).
"""
abstract type AbstractDensityResponse end

# --- Concrete types ---

"""
    RickerDensity(α, β)

Ricker-type overcompensatory density dependence:
    f(N) = exp(α - β × N)

At N=0, f = exp(α). Decays exponentially with density.
Can produce overcompensation (f < 1 at high densities even if α > 0).
"""
struct RickerDensity{T<:Real} <: AbstractDensityResponse
    α::T
    β::T
end

RickerDensity(; α=0.0, β=1.0) = RickerDensity(α, β)

(f::RickerDensity)(N::Real) = exp(f.α - f.β * N)

"""
    BevertonHoltDensity(α, β)

Beverton-Holt compensatory density dependence:
    f(N) = α / (1 + β × N)

At N=0, f = α. Monotonically decreasing, bounded below by 0.
Classic contest competition form — carrying capacity K = (α-1)/β when α > 1.
"""
struct BevertonHoltDensity{T<:Real} <: AbstractDensityResponse
    α::T
    β::T
end

BevertonHoltDensity(; α=1.0, β=1.0) = BevertonHoltDensity(α, β)

(f::BevertonHoltDensity)(N::Real) = f.α / (1 + f.β * N)

"""
    UsherDensity(α, β)

Usher/logistic-type sigmoidal density dependence:
    f(N) = 1 / (1 + exp(α + β × N))

Logistic (sigmoidal) decline with density. At low N, f ≈ 1/(1+exp(α)).
"""
struct UsherDensity{T<:Real} <: AbstractDensityResponse
    α::T
    β::T
end

UsherDensity(; α=0.0, β=1.0) = UsherDensity(α, β)

(f::UsherDensity)(N::Real) = 1.0 / (1.0 + exp(f.α + f.β * N))

"""
    LogisticDensity(r, K)

Linear logistic density dependence:
    f(N) = max(0, 1 - N/K)

At N=0, f = 1. Linearly decreasing to 0 at N = K.
The parameter `r` is unused in the modifier but conventionally stored
for the full logistic growth context.
"""
struct LogisticDensity{T<:Real} <: AbstractDensityResponse
    r::T
    K::T
end

LogisticDensity(; r=1.0, K=100.0) = LogisticDensity(r, K)

(f::LogisticDensity)(N::Real) = max(zero(N), 1 - N / f.K)

"""
    ThetaLogisticDensity(r, K, θ)

Theta-logistic (generalized) density dependence:
    f(N) = max(0, 1 - (N/K)^θ)

Generalizes linear logistic (θ=1). θ < 1 gives concave decline (strong
regulation at low N), θ > 1 gives convex decline (weak regulation until
near K).
"""
struct ThetaLogisticDensity{T<:Real} <: AbstractDensityResponse
    r::T
    K::T
    θ::T
end

ThetaLogisticDensity(; r=1.0, K=100.0, θ=1.0) = ThetaLogisticDensity(r, K, θ)

(f::ThetaLogisticDensity)(N::Real) = max(zero(N), 1 - (N / f.K)^f.θ)

"""
    ConstantDensity()

No density dependence (identity modifier). Always returns 1.0.
Useful as a placeholder when some vital rates are density-independent.
"""
struct ConstantDensity <: AbstractDensityResponse end

(::ConstantDensity)(N::Real) = one(N)

# --- Density specification for MPMs ---

"""
    DensityVitalRateSpec

Specifies density-dependent modification of vital rates in a projection matrix.
Maps each vital rate process to a density response function with optional time delay.

# Fields
- `survival`: Density response for survival rates
- `growth`: Density response for growth/transition rates
- `fecundity`: Density response for fecundity rates
- `recruitment`: Density response for recruitment rates
- `time_delay`: Number of time steps to lag the density signal (default 1 = current)

# Example
```julia
spec = DensityVitalRateSpec(
    survival = BevertonHoltDensity(α=1.0, β=0.01),
    fecundity = RickerDensity(α=0.0, β=0.005),
)
```
"""
struct DensityVitalRateSpec{S<:AbstractDensityResponse,
                           G<:AbstractDensityResponse,
                           F<:AbstractDensityResponse,
                           R<:AbstractDensityResponse}
    survival::S
    growth::G
    fecundity::F
    recruitment::R
    time_delay::Int
end

function DensityVitalRateSpec(;
        survival::AbstractDensityResponse = ConstantDensity(),
        growth::AbstractDensityResponse = ConstantDensity(),
        fecundity::AbstractDensityResponse = ConstantDensity(),
        recruitment::AbstractDensityResponse = ConstantDensity(),
        time_delay::Int = 1)
    DensityVitalRateSpec(survival, growth, fecundity, recruitment, time_delay)
end

"""
    apply_density(spec::DensityVitalRateSpec, matU::AbstractMatrix,
                  matF::AbstractMatrix, N::Real)

Apply density-dependent modification to survival (U) and fecundity (F) matrices.
`growth`, `recruitment`, and lagged density responses are not yet implemented for
MPM matrix transforms and will raise an error if requested.
Returns modified (U_dd, F_dd) matrices.
"""
function _validate_density_spec(spec::DensityVitalRateSpec)
    spec.growth isa ConstantDensity ||
        throw(ArgumentError("DensityVitalRateSpec.growth is not applied by apply_density; use ConstantDensity() or transform the transition matrix explicitly."))
    spec.recruitment isa ConstantDensity ||
        throw(ArgumentError("DensityVitalRateSpec.recruitment is not applied by apply_density; use ConstantDensity() or transform the recruitment matrix explicitly."))
    spec.time_delay == 1 ||
        throw(ArgumentError("DensityVitalRateSpec.time_delay is not applied by apply_density; only time_delay = 1 is currently supported."))
    return spec
end

function apply_density(spec::DensityVitalRateSpec,
                       matU::AbstractMatrix, matF::AbstractMatrix, N::Real)
    _validate_density_spec(spec)
    surv_mod = spec.survival(N)
    fec_mod = spec.fecundity(N)
    return matU .* surv_mod, matF .* fec_mod
end

"""
    apply_density(spec::DensityVitalRateSpec, matA::AbstractMatrix, N::Real;
                  fec_rows=1:1)

Apply density dependence to a combined A matrix. Elements in `fec_rows`
are modified by the fecundity response; others by the survival response.
`growth`, `recruitment`, and lagged density responses are not yet implemented for
combined matrices and will raise an error if requested.
"""
function apply_density(spec::DensityVitalRateSpec,
                       matA::AbstractMatrix, N::Real;
                       fec_rows::AbstractVector{Int}=1:1)
    _validate_density_spec(spec)
    result = copy(float.(matA))
    surv_mod = spec.survival(N)
    fec_mod = spec.fecundity(N)

    for j in axes(result, 2), i in axes(result, 1)
        if i in fec_rows
            result[i, j] *= fec_mod
        else
            result[i, j] *= surv_mod
        end
    end
    return result
end
