"""
Generation time calculations.
Based on Rage::gen_time with 3 methods.
"""

"""
    gen_time(matU::AbstractMatrix, matR::AbstractMatrix;
             start::Int=1, method::Symbol=:R0)

Generation time: mean age of parents of offspring produced by a cohort.

# Methods
- `:R0`: T = log(R0) / log(λ) — most common definition
- `:cohort`: T = sum(x * lx * mx) / R0 — cohort-based mean age of reproduction
- `:age_diff`: T = mean age of parents - mean age of offspring (from stable population)

# Arguments
- `matU`: Survival matrix
- `matR`: Reproductive matrix
- `start`: Starting stage
"""
function gen_time(matU::AbstractMatrix, matR::AbstractMatrix;
                  start::Int=1, method::Symbol=:R0)
    n = size(matU, 1)
    _validate_start(start, n)

    if method == :R0
        R0 = net_repro_rate(matU, matR; start=start)
        R0 <= 0 && return Inf
        A = matU .+ matR
        λ = lambda(A)
        λ <= 0 && return Inf
        return log(R0) / log(λ)

    elseif method == :cohort
        lx = mpm_to_lx(matU; start=start)
        mx = mpm_to_mx(matU, matR; start=start)
        nages = min(length(lx), length(mx))
        R0 = sum(lx[i] * mx[i] for i in 1:nages)
        R0 <= 0 && return Inf
        # Mean age of reproduction
        return sum((i - 1) * lx[i] * mx[i] for i in 1:nages) / R0

    elseif method == :age_diff
        A = matU .+ matR
        λ = lambda(A)
        λ <= 0 && return Inf
        lx = mpm_to_lx(matU; start=start)
        mx = mpm_to_mx(matU, matR; start=start)
        nages = min(length(lx), length(mx))
        # T = Σ x * lx * mx * λ^(-x) / Σ lx * mx * λ^(-x)
        num = sum((i - 1) * lx[i] * mx[i] * λ^(-(i - 1)) for i in 1:nages)
        den = sum(lx[i] * mx[i] * λ^(-(i - 1)) for i in 1:nages)
        den ≈ 0 && return Inf
        return num / den

    else
        throw(ArgumentError("method must be :R0, :cohort, or :age_diff"))
    end
end
