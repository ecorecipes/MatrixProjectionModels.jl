"""
Shape measures for survivorship and reproduction.
Based on Rage::shape_surv, shape_rep.
"""

"""
    shape_surv(lx::AbstractVector)

Shape of survivorship curve.
Ranges from -1/6 (Type III, high early mortality) to +1/6 (Type I, late mortality).
0 indicates constant mortality (Type II).

Calculated as the area under the scaled survivorship curve minus 0.5:
shape = ∫₀¹ l(x/ω) dx - 0.5, scaled to [-1/6, 1/6].
"""
function shape_surv(lx::AbstractVector)
    n = length(lx)
    n < 2 && return 0.0

    # Normalize lx to start at 1 and have unit x-range
    lx_norm = lx ./ lx[1]
    x = range(0, 1; length=n)

    # Area under normalized curve using trapezoidal rule
    area = area_under_curve(collect(x), lx_norm)

    # Shape: deviation from type II (area = 0.5)
    return area - 0.5
end

"""
    shape_rep(matU::AbstractMatrix, matR::AbstractMatrix;
              start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)

Shape of reproduction curve (analogous to shape_surv).
Positive values indicate reproduction concentrated late in life.
Negative values indicate early reproduction.
"""
function shape_rep(matU::AbstractMatrix, matR::AbstractMatrix;
                   start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)
    mx = mpm_to_mx(matU, matR; start=start, xmax=xmax, lx_crit=lx_crit)

    n = length(mx)
    n < 2 && return 0.0

    # Normalize to cumulative reproduction (like survivorship but for reproduction)
    cum_mx = cumsum(mx)
    total = cum_mx[end]
    total ≈ 0 && return 0.0
    cum_mx ./= total

    x = range(0, 1; length=n)
    area = area_under_curve(collect(x), cum_mx)

    return area - 0.5
end
