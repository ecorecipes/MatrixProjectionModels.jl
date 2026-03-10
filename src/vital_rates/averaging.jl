"""
Weighted averaging of vital rates across stages.
Based on Rage::vr_survival, vr_growth, vr_fecundity, etc.
"""

"""
    vr_survival(matU::AbstractMatrix; weights::Union{Nothing,AbstractVector}=nothing)

Mean survival rate across stages, optionally weighted.
Default weights: stable stage distribution from U.
"""
function vr_survival(matU::AbstractMatrix; weights::Union{Nothing,AbstractVector}=nothing)
    vals = vr_vec_survival(matU)
    w = _get_weights(matU, weights)
    return dot(vals, w) / sum(w)
end

"""
    vr_growth(matU::AbstractMatrix; weights::Union{Nothing,AbstractVector}=nothing)

Mean growth (progression) rate across stages.
"""
function vr_growth(matU::AbstractMatrix; weights::Union{Nothing,AbstractVector}=nothing)
    vals = vr_vec_growth(matU)
    w = _get_weights(matU, weights)
    return dot(vals, w) / sum(w)
end

"""
    vr_shrinkage(matU::AbstractMatrix; weights::Union{Nothing,AbstractVector}=nothing)

Mean shrinkage (retrogression) rate across stages.
"""
function vr_shrinkage(matU::AbstractMatrix; weights::Union{Nothing,AbstractVector}=nothing)
    vals = vr_vec_shrinkage(matU)
    w = _get_weights(matU, weights)
    return dot(vals, w) / sum(w)
end

"""
    vr_stasis(matU::AbstractMatrix; weights::Union{Nothing,AbstractVector}=nothing)

Mean stasis rate across stages.
"""
function vr_stasis(matU::AbstractMatrix; weights::Union{Nothing,AbstractVector}=nothing)
    vals = vr_vec_stasis(matU)
    w = _get_weights(matU, weights)
    return dot(vals, w) / sum(w)
end

"""
    vr_fecundity(matU::AbstractMatrix, matF::AbstractMatrix;
                 weights::Union{Nothing,AbstractVector}=nothing)

Mean fecundity rate across stages, conditional on survival.
"""
function vr_fecundity(matU::AbstractMatrix, matF::AbstractMatrix;
                      weights::Union{Nothing,AbstractVector}=nothing)
    vals = vr_vec_reproduction(matU, matF)
    w = _get_weights(matU, weights)
    return dot(vals, w) / sum(w)
end

# Helper: default weights from stable distribution of U
function _get_weights(matU::AbstractMatrix, weights::Union{Nothing,AbstractVector})
    if weights !== nothing
        return weights
    end
    # Use uniform weights if U is all zeros
    n = size(matU, 1)
    col_sums = vec(sum(matU; dims=1))
    if all(col_sums .== 0)
        return ones(n) / n
    end
    # Use stable distribution of A = U (ignoring reproduction for weighting)
    w = stable_distribution(matU)
    # Ensure non-negative
    w = max.(w, 0.0)
    s = sum(w)
    return s > 0 ? w : ones(n) / n
end
