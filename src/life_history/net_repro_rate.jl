"""
Net reproductive rate (R0).
Based on Rage::net_repro_rate.
"""

"""
    net_repro_rate(matU::AbstractMatrix, matR::AbstractMatrix;
                   start::Int=1, method::Symbol=:generation)

Net reproductive rate R0: expected number of offspring per individual
over its lifetime.

# Methods
- `:generation`: R0 = sum(lx .* mx) from life table
- `:fundamental`: R0 from fundamental matrix: R0 = (R * N)[start, start]

# Arguments
- `matU`: Survival/growth transition matrix
- `matR`: Reproductive matrix (F, or F+C)
- `start`: Starting stage
- `method`: Computation method
"""
function net_repro_rate(matU::AbstractMatrix, matR::AbstractMatrix;
                        start::Int=1, method::Symbol=:generation)
    n = size(matU, 1)
    _validate_start(start, n)

    if method == :generation
        lx = mpm_to_lx(matU; start=start)
        mx = mpm_to_mx(matU, matR; start=start)
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
