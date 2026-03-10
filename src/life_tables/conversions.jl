"""
Life table schedule conversions.
Based on Rage life table conversion utilities.
"""

"""
    lx_to_px(lx::AbstractVector)

Convert survivorship (lx) to survival probability (px).
px[x] = lx[x+1] / lx[x]. Last value is 0 (absorbing).
"""
function lx_to_px(lx::AbstractVector)
    n = length(lx)
    px = zeros(n)
    for i in 1:(n-1)
        px[i] = lx[i] > 0 ? lx[i+1] / lx[i] : 0.0
    end
    px[n] = 0.0
    return px
end

"""
    lx_to_hx(lx::AbstractVector)

Convert survivorship (lx) to hazard rate (hx).
hx[x] = -log(px[x]). Last value is Inf (absorbing).
"""
function lx_to_hx(lx::AbstractVector)
    px = lx_to_px(lx)
    return px_to_hx(px)
end

"""
    px_to_lx(px::AbstractVector)

Convert survival probability (px) to survivorship (lx).
lx[1] = 1, lx[x+1] = lx[x] * px[x].
"""
function px_to_lx(px::AbstractVector)
    n = length(px)
    lx = ones(n + 1)
    for i in 1:n
        lx[i+1] = lx[i] * px[i]
    end
    return lx
end

"""
    px_to_hx(px::AbstractVector)

Convert survival probability (px) to hazard rate (hx).
hx[x] = -log(px[x]).
"""
function px_to_hx(px::AbstractVector)
    hx = similar(px)
    for i in eachindex(px)
        hx[i] = px[i] > 0 ? -log(px[i]) : Inf
    end
    return hx
end

"""
    hx_to_lx(hx::AbstractVector)

Convert hazard rate (hx) to survivorship (lx).
lx[1] = 1, lx[x+1] = lx[x] * exp(-hx[x]).
"""
function hx_to_lx(hx::AbstractVector)
    n = length(hx)
    lx = ones(n + 1)
    for i in 1:n
        lx[i+1] = lx[i] * exp(-hx[i])
    end
    return lx
end

"""
    hx_to_px(hx::AbstractVector)

Convert hazard rate (hx) to survival probability (px).
px[x] = exp(-hx[x]).
"""
function hx_to_px(hx::AbstractVector)
    return exp.(-hx)
end
