"""
Longevity calculation.
Based on Rage::longevity.
"""

"""
    longevity(matU::AbstractMatrix; start::Int=1, lx_crit::Real=0.01)

Age at which survivorship lx first drops below `lx_crit`.
This is an estimate of maximum lifespan.
"""
function longevity(matU::AbstractMatrix; start::Int=1, lx_crit::Real=0.01, xmax::Int=1000)
    n = size(matU, 1)
    _validate_start(start, n)
    n_vec = zeros(Float64, n)
    n_vec[start] = 1.0
    tmp = similar(n_vec)
    for x in 1:xmax
        mul!(tmp, matU, n_vec)
        n_vec, tmp = tmp, n_vec
        sum(n_vec) < lx_crit && return x
    end
    return xmax
end
