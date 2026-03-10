"""
Mortality models as callable structs returning hazard rate h(x).
Based on mpmsim R package.
"""

abstract type AbstractMortalityModel end

"""
    GompertzMortality(b0, b1)

Gompertz mortality: h(x) = b0 * exp(b1 * x)
"""
struct GompertzMortality{T<:Real} <: AbstractMortalityModel
    b0::T
    b1::T
end
GompertzMortality(b0, b1) = GompertzMortality(promote(b0, b1)...)
(m::GompertzMortality)(x) = m.b0 * exp(m.b1 * x)

"""
    GompertzMakehamMortality(b0, b1, C)

Gompertz-Makeham mortality: h(x) = b0 * exp(b1 * x) + C
"""
struct GompertzMakehamMortality{T<:Real} <: AbstractMortalityModel
    b0::T
    b1::T
    C::T
end
GompertzMakehamMortality(b0, b1, C) = GompertzMakehamMortality(promote(b0, b1, C)...)
(m::GompertzMakehamMortality)(x) = m.b0 * exp(m.b1 * x) + m.C

"""
    ExponentialMortality(C)

Exponential (constant) mortality: h(x) = C
"""
struct ExponentialMortality{T<:Real} <: AbstractMortalityModel
    C::T
end
(m::ExponentialMortality)(x) = m.C

"""
    SilerMortality(a0, a1, C, b0, b1)

Siler mortality: h(x) = a0 * exp(-a1 * x) + C + b0 * exp(b1 * x)
Bathtub-shaped: juvenile + constant + senescent components.
"""
struct SilerMortality{T<:Real} <: AbstractMortalityModel
    a0::T
    a1::T
    C::T
    b0::T
    b1::T
end
SilerMortality(a0, a1, C, b0, b1) = SilerMortality(promote(a0, a1, C, b0, b1)...)
(m::SilerMortality)(x) = m.a0 * exp(-m.a1 * x) + m.C + m.b0 * exp(m.b1 * x)

"""
    WeibullMortality(b0, b1)

Weibull mortality: h(x) = b0 * b1 * (b1 * x)^(b0 - 1)
"""
struct WeibullMortality{T<:Real} <: AbstractMortalityModel
    b0::T
    b1::T
end
WeibullMortality(b0, b1) = WeibullMortality(promote(b0, b1)...)
(m::WeibullMortality)(x) = m.b0 * m.b1 * (m.b1 * x)^(m.b0 - 1)

"""
    WeibullMakehamMortality(b0, b1, C)

Weibull-Makeham mortality: h(x) = b0 * b1 * (b1 * x)^(b0 - 1) + C
"""
struct WeibullMakehamMortality{T<:Real} <: AbstractMortalityModel
    b0::T
    b1::T
    C::T
end
WeibullMakehamMortality(b0, b1, C) = WeibullMakehamMortality(promote(b0, b1, C)...)
(m::WeibullMakehamMortality)(x) = m.b0 * m.b1 * (m.b1 * x)^(m.b0 - 1) + m.C

"""
    model_survival(model::AbstractMortalityModel; ages=0:1000, truncate=0.01)

Compute survival schedule from a mortality model.
Returns NamedTuple (x, hx, lx, qx, px) where:
- x: ages
- hx: hazard rate
- lx: survivorship (cumulative survival)
- qx: mortality probability per interval
- px: survival probability per interval

Truncates at age where lx < `truncate`.
"""
function model_survival(model::AbstractMortalityModel;
                        ages::AbstractVector=0:1000, truncate::Real=0.01)
    x = collect(ages)
    hx = [model(Float64(a)) for a in x]

    # Cumulative hazard → survivorship
    # lx = exp(-∫₀ˣ h(t) dt), approximate with cumulative sum
    Δ = length(x) > 1 ? Float64(x[2] - x[1]) : 1.0
    cum_h = cumsum(hx) .* Δ
    lx = exp.(-cum_h)

    # Truncate
    last_idx = findfirst(l -> l < truncate, lx)
    if last_idx !== nothing
        last_idx = max(last_idx, 2)  # Keep at least 2 points
        x = x[1:last_idx]
        hx = hx[1:last_idx]
        lx = lx[1:last_idx]
    end

    # Derived quantities
    n = length(x)
    px = zeros(n)
    qx = zeros(n)
    for i in 1:(n-1)
        px[i] = lx[i] > 0 ? lx[i+1] / lx[i] : 0.0
        qx[i] = 1.0 - px[i]
    end
    # Last age class: absorbing
    px[n] = 0.0
    qx[n] = 1.0

    return (x=x, hx=hx, lx=lx, qx=qx, px=px)
end
