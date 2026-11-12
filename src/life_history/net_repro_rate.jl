"""
Net reproductive rate (R0).
Based on Rage::net_repro_rate.
"""

"""
    net_repro_rate(matU::AbstractMatrix, matR::AbstractMatrix;
                   start::Int=1, method::Symbol=:generation,
                   xmax::Int=1000, lx_crit::Real=0.01)

Net reproductive rate R0: expected number of offspring per individual
over its lifetime.

# Methods
- `:generation`: R0 = sum(lx .* mx) from the life table. The life table is truncated at
  the first age at which survivorship falls below `lx_crit` (or at `xmax`), so the result
  is a lower bound that omits reproduction by the last `lx_crit` of a cohort. For a life
  cycle with a long-lived reproductive stage the shortfall can be a few percent with the
  default `lx_crit`; lower `lx_crit` to converge on the `:fundamental` value.
- `:fundamental`: R0 from the fundamental matrix, R0 = (R * N)[start, start]. Exact, with
  no truncation; `xmax` and `lx_crit` are ignored.

# Arguments
- `matU`: Survival/growth transition matrix
- `matR`: Reproductive matrix (F, or F+C)
- `start`: Starting stage
- `method`: Computation method
- `xmax`: Maximum age of the life table (`:generation` only)
- `lx_crit`: Survivorship below which the life table stops (`:generation` only)
"""
function net_repro_rate(matU::AbstractMatrix, matR::AbstractMatrix;
                        start::Int=1, method::Symbol=:generation,
                        xmax::Int=1000, lx_crit::Real=0.01)
    n = size(matU, 1)
    _validate_start(start, n)

    if method == :generation
        lx = mpm_to_lx(matU; start=start, xmax=xmax, lx_crit=lx_crit)
        mx = mpm_to_mx(matU, matR; start=start, xmax=xmax, lx_crit=lx_crit)
        # Ensure same length
        nages = min(length(lx), length(mx))
        return sum(lx[i] * mx[i] for i in 1:nages)
    elseif method == :fundamental
        N = _fundamental_matrix(matU)
        N === nothing && return Inf
        R0_mat = matR * N
        return R0_mat[start, start]
    else
        throw(ArgumentError("method must be :generation or :fundamental"))
    end
end
