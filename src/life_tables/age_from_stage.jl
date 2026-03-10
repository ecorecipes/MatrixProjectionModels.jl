"""
Age-from-stage life table calculations.
Based on Rage::mpm_to_lx, mpm_to_px, mpm_to_hx, mpm_to_mx, mpm_to_table.

These functions convert stage-classified MPMs to age-classified demographic
schedules by repeated matrix multiplication.
"""

"""
    mpm_to_lx(matU::AbstractMatrix; start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)

Calculate age-specific survivorship (lx) from a survival matrix U.
lx[x] is the probability of surviving from the start stage to age x.

# Arguments
- `matU`: Survival/growth transition matrix (n × n)
- `start`: Index of the starting stage (default: 1)
- `xmax`: Maximum age to compute (default: 1000)
- `lx_crit`: Stop when lx drops below this value (default: 0.01)

# Returns
Vector of survivorship values, starting with lx[1] = 1.0.
"""
function mpm_to_lx(matU::AbstractMatrix; start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)
    n = size(matU, 1)
    _validate_start(start, n)

    lx = Vector{Float64}(undef, xmax + 1)
    lx[1] = 1.0
    # State vector: unit vector at start stage
    n_vec = zeros(Float64, n)
    n_vec[start] = 1.0
    tmp = similar(n_vec)
    count = 1

    for x in 1:xmax
        mul!(tmp, matU, n_vec)
        n_vec, tmp = tmp, n_vec
        surv = sum(n_vec)
        count += 1
        lx[count] = surv
        surv < lx_crit && break
    end

    return resize!(lx, count)
end

"""
    mpm_to_px(matU::AbstractMatrix; start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)

Age-specific survival probability: px[x] = lx[x+1] / lx[x].
"""
function mpm_to_px(matU::AbstractMatrix; start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)
    lx = mpm_to_lx(matU; start=start, xmax=xmax, lx_crit=lx_crit)
    return lx_to_px(lx)
end

"""
    mpm_to_hx(matU::AbstractMatrix; start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)

Age-specific hazard rate: hx[x] = -log(px[x]).
"""
function mpm_to_hx(matU::AbstractMatrix; start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)
    lx = mpm_to_lx(matU; start=start, xmax=xmax, lx_crit=lx_crit)
    return lx_to_hx(lx)
end

"""
    mpm_to_mx(matU::AbstractMatrix, matR::AbstractMatrix;
              start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)

Age-specific reproduction (mx).
matR is the reproductive matrix (F + C, or just F).

mx[x] = expected reproduction at age x, conditional on being alive.
"""
function mpm_to_mx(matU::AbstractMatrix, matR::AbstractMatrix;
                   start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)
    n = size(matU, 1)
    _validate_start(start, n)

    lx = mpm_to_lx(matU; start=start, xmax=xmax, lx_crit=lx_crit)
    nages = length(lx)

    mx = zeros(nages)
    n_vec = zeros(Float64, n)
    n_vec[start] = 1.0

    for x in 1:nages
        # Reproduction at age x (conditional on surviving to x)
        r_vec = matR * n_vec
        mx[x] = lx[x] > 0 ? sum(r_vec) / lx[x] : 0.0
        # Advance age
        n_vec = matU * n_vec
    end

    return mx
end

"""
    mpm_to_table(matU::AbstractMatrix, matR::AbstractMatrix;
                 start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)

Compute a full life table from stage-classified matrices.

# Returns
NamedTuple (x, lx, px, hx, mx) where:
- x: age (0-indexed)
- lx: survivorship
- px: survival probability
- hx: hazard rate
- mx: age-specific reproduction
"""
function mpm_to_table(matU::AbstractMatrix, matR::AbstractMatrix;
                      start::Int=1, xmax::Int=1000, lx_crit::Real=0.01)
    lx = mpm_to_lx(matU; start=start, xmax=xmax, lx_crit=lx_crit)
    px = lx_to_px(lx)
    hx = lx_to_hx(lx)
    mx = mpm_to_mx(matU, matR; start=start, xmax=xmax, lx_crit=lx_crit)

    # Align lengths (mx may differ slightly)
    nages = length(lx)
    if length(mx) < nages
        append!(mx, zeros(nages - length(mx)))
    elseif length(mx) > nages
        mx = mx[1:nages]
    end

    x = collect(0:(nages - 1))
    return (x=x, lx=lx, px=px, hx=hx, mx=mx)
end
