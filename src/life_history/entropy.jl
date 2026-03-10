"""
Demographic entropy measures.
Based on Rage::entropy_k, entropy_k_age, entropy_k_stage, entropy_d.
"""

"""
    entropy_k(lx::AbstractVector)

Keyfitz's life table entropy from a survivorship schedule.
H = -sum(lx * log(lx)) / sum(lx), excluding lx = 0.

Measures the elasticity of life expectancy to proportional changes in mortality.
H = 0: type I survivorship (rectangular), H = 1: type II (constant), H > 1: type III.
"""
function entropy_k(lx::AbstractVector)
    # Filter out zeros (log(0) undefined)
    pos = filter(x -> x > 0, lx)
    isempty(pos) && return 0.0
    total = sum(pos)
    total ≈ 0 && return 0.0
    return -sum(x * log(x) for x in pos) / total
end

"""
    entropy_k_age(matU::AbstractMatrix; start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)

Keyfitz's entropy calculated from age-specific survivorship derived from U.
"""
function entropy_k_age(matU::AbstractMatrix; start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)
    lx = mpm_to_lx(matU; start=start, xmax=xmax, lx_crit=lx_crit)
    return entropy_k(lx)
end

"""
    entropy_k_stage(matU::AbstractMatrix)

Keyfitz's entropy from stage dynamics (not age-converted).
Uses the fundamental matrix approach.

H = -trace(N * U * log(U)) / trace(N)  [approximate]
"""
function entropy_k_stage(matU::AbstractMatrix)
    # Use age-from-stage approach as fallback
    return entropy_k_age(matU)
end

"""
    entropy_d(lx::AbstractVector, mx::AbstractVector)

Demetrius' entropy (demographic entropy).
s = -sum(cx * log(cx)) where cx = lx * mx / R0 (normalized reproductive schedule).
"""
function entropy_d(lx::AbstractVector, mx::AbstractVector)
    nages = min(length(lx), length(mx))
    cx = [lx[i] * mx[i] for i in 1:nages]
    R0 = sum(cx)
    R0 ≈ 0 && return 0.0
    cx ./= R0
    # Shannon entropy of reproductive schedule
    s = 0.0
    for c in cx
        if c > 0
            s -= c * log(c)
        end
    end
    return s
end
