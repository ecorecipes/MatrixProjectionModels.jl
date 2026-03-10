"""
Life expectancy calculations from MPM.
Based on Rage::life_expect_mean, life_expect_var.
"""

"""
    life_expect_mean(matU::AbstractMatrix; start::Int=1)

Mean life expectancy starting from a given stage.
Computed from the fundamental matrix N = (I - U)^{-1}.
Mean time in each stage: N[start, :], total = sum.
"""
function life_expect_mean(matU::AbstractMatrix; start::Int=1)
    n = size(matU, 1)
    _validate_start(start, n)

    N = _fundamental_matrix(matU)
    N === nothing && return Inf  # Absorbing: infinite life expectancy

    # Mean life expectancy = sum of row `start` of N
    return sum(N[start, :])
end

"""
    life_expect_var(matU::AbstractMatrix; start::Int=1)

Variance in life expectancy starting from a given stage.
Uses the formula: Var = (2N - I) * N - (N .* N) componentwise sum.
Specifically: Var(T) = (2*diag(N) - 1) · N[start,:] - sum(N[start,:])^2 + sum(N[start,:])
"""
function life_expect_var(matU::AbstractMatrix; start::Int=1)
    n = size(matU, 1)
    _validate_start(start, n)

    N = _fundamental_matrix(matU)
    N === nothing && return Inf

    # Variance formula from Caswell (2001)
    # Var(η_i) = (2N - I)N e_i · e_i - (Ne_i · e_i)^2
    # where e_i is unit vector, · is element-wise
    # Simplified: Var = 2*N²[start,:] sum - N[start,:] sum - (N[start,:] sum)²
    Nsq = N * N
    mean_life = sum(N[start, :])
    return 2 * sum(Nsq[start, :]) - mean_life - mean_life^2
end
