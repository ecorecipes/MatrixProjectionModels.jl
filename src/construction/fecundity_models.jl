"""
Fecundity models as callable structs returning f(x).
Based on mpmsim R package.
"""

abstract type AbstractFecundityModel end

"""
    LogisticFecundity(A, k, x_mid)

Logistic fecundity: f(x) = A / (1 + exp(-k * (x - x_mid)))
"""
struct LogisticFecundity{T<:Real} <: AbstractFecundityModel
    A::T       # Maximum fecundity
    k::T       # Steepness
    x_mid::T   # Midpoint age
end
LogisticFecundity(A, k, x_mid) = LogisticFecundity(promote(A, k, x_mid)...)
(m::LogisticFecundity)(x) = m.A / (1 + exp(-m.k * (x - m.x_mid)))

"""
    StepFecundity(A, x_onset)

Step fecundity: f(x) = A if x >= x_onset, else 0.
"""
struct StepFecundity{T<:Real} <: AbstractFecundityModel
    A::T       # Fecundity level
    x_onset::T # Onset age
end
StepFecundity(A, x_onset) = StepFecundity(promote(A, x_onset)...)
(m::StepFecundity)(x) = x >= m.x_onset ? m.A : zero(m.A)

"""
    VonBertalanffyFecundity(A, k)

Von Bertalanffy fecundity: f(x) = A * (1 - exp(-k * x))
"""
struct VonBertalanffyFecundity{T<:Real} <: AbstractFecundityModel
    A::T    # Asymptotic fecundity
    k::T    # Growth rate
end
VonBertalanffyFecundity(A, k) = VonBertalanffyFecundity(promote(A, k)...)
(m::VonBertalanffyFecundity)(x) = m.A * (1 - exp(-m.k * x))

"""
    NormalFecundity(A, mu, sd)

Normal (Gaussian) fecundity: f(x) = A * exp(-((x - mu)^2) / (2 * sd^2))
"""
struct NormalFecundity{T<:Real} <: AbstractFecundityModel
    A::T    # Peak fecundity
    mu::T   # Peak age
    sd::T   # Spread
end
NormalFecundity(A, mu, sd) = NormalFecundity(promote(A, mu, sd)...)
(m::NormalFecundity)(x) = m.A * exp(-((x - m.mu)^2) / (2 * m.sd^2))

"""
    HadwigerFecundity(a, b, C)

Hadwiger fecundity: f(x) = (a * b / C) * (C/x)^(3/2) * exp(-b^2 * (C/x + x/C - 2))
A model of human fertility.
"""
struct HadwigerFecundity{T<:Real} <: AbstractFecundityModel
    a::T
    b::T
    C::T
end
HadwigerFecundity(a, b, C) = HadwigerFecundity(promote(a, b, C)...)
function (m::HadwigerFecundity)(x)
    x <= 0 && return zero(m.a)
    r = m.C / x
    return (m.a * m.b / m.C) * r^(3/2) * exp(-m.b^2 * (r + x / m.C - 2))
end

"""
    model_fecundity(model::AbstractFecundityModel; ages=0:100)

Compute fecundity schedule from a fecundity model.
Returns NamedTuple (x, fx).
"""
function model_fecundity(model::AbstractFecundityModel; ages::AbstractVector=0:100)
    x = collect(ages)
    fx = [model(Float64(a)) for a in x]
    return (x=x, fx=fx)
end
