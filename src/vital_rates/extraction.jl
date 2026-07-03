"""
Vital rate extraction from MPM matrices.
Based on Rage::vr_vec_* functions.
"""

"""
    vr_vec_survival(matU::AbstractMatrix)

Stage-specific survival: column sums of U.
"""
function vr_vec_survival(matU::AbstractMatrix)
    return vec(sum(matU; dims=1))
end

"""
    vr_vec_growth(matU::AbstractMatrix)

Stage-specific growth (progression) probability.
For each stage j, the probability of transitioning to a higher stage,
conditional on surviving.
"""
function vr_vec_growth(matU::AbstractMatrix)
    n = size(matU, 1)
    growth = zeros(n)
    for j in 1:n
        col_sum = sum(matU[:, j])
        if col_sum > 0
            # Sum of transitions to higher stages
            higher = sum(matU[i, j] for i in (j+1):n; init=0.0)
            growth[j] = higher / col_sum
        end
    end
    return growth
end

"""
    vr_vec_shrinkage(matU::AbstractMatrix)

Stage-specific shrinkage (retrogression) probability.
For each stage j, the probability of transitioning to a lower stage,
conditional on surviving.
"""
function vr_vec_shrinkage(matU::AbstractMatrix)
    n = size(matU, 1)
    shrinkage = zeros(n)
    for j in 1:n
        col_sum = sum(matU[:, j])
        if col_sum > 0
            # Sum of transitions to lower stages
            lower = sum(matU[i, j] for i in 1:(j-1); init=0.0)
            shrinkage[j] = lower / col_sum
        end
    end
    return shrinkage
end

"""
    vr_vec_stasis(matU::AbstractMatrix)

Stage-specific stasis probability.
For each stage j, the probability of remaining in the same stage,
conditional on surviving.
"""
function vr_vec_stasis(matU::AbstractMatrix)
    n = size(matU, 1)
    stasis = zeros(n)
    for j in 1:n
        col_sum = sum(matU[:, j])
        if col_sum > 0
            stasis[j] = matU[j, j] / col_sum
        end
    end
    return stasis
end

"""
    vr_vec_reproduction(matU::AbstractMatrix, matR::AbstractMatrix)

Stage-specific reproduction rate, conditional on surviving.
matR is the reproductive matrix (F, or F + C).
"""
function vr_vec_reproduction(matU::AbstractMatrix, matR::AbstractMatrix)
    n = size(matU, 1)
    repro = zeros(n)
    surv = vr_vec_survival(matU)
    for j in 1:n
        r = sum(matR[:, j])
        repro[j] = surv[j] > 0 ? r / surv[j] : 0.0
    end
    return repro
end

"""
    vr_vec_dorm_enter(matU::AbstractMatrix; dorm_stages::AbstractVector{Int})

Stage-specific probability of entering dormancy.
"""
function vr_vec_dorm_enter(matU::AbstractMatrix; dorm_stages::AbstractVector{Int})
    n = size(matU, 1)
    active = setdiff(1:n, dorm_stages)
    result = zeros(n)
    for j in active
        col_sum = sum(matU[:, j])
        if col_sum > 0
            dorm_total = sum(matU[i, j] for i in dorm_stages; init=0.0)
            result[j] = dorm_total / col_sum
        end
    end
    return result
end

"""
    vr_vec_dorm_exit(matU::AbstractMatrix; dorm_stages::AbstractVector{Int})

Stage-specific probability of exiting dormancy.
"""
function vr_vec_dorm_exit(matU::AbstractMatrix; dorm_stages::AbstractVector{Int})
    n = size(matU, 1)
    active = setdiff(1:n, dorm_stages)
    result = zeros(n)
    for j in dorm_stages
        col_sum = sum(matU[:, j])
        if col_sum > 0
            active_total = sum(matU[i, j] for i in active; init=0.0)
            result[j] = active_total / col_sum
        end
    end
    return result
end
