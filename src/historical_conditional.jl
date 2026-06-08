"""
Historical conditional extraction from lagged (historical) MPMs.

Extracts ahistorical conditional matrices from a historical (memory/lag)
matrix, conditioned on the previous state.

Reference: lefko3's cond_hmpm() function (Shefferson & Ehrlen).
"""

"""
    conditional_ahistorical(lag_matrix::AbstractMatrix, n_stages::Int;
                            condition_stage::Union{Int,Nothing}=nothing)

Extract ahistorical conditional matrices from a historical (expanded lag) matrix.

A historical MPM of size (n²×n²) encodes transitions (stage_t-1, stage_t) → (stage_t, stage_t+1).
This function extracts the n×n ahistorical matrix for a specific previous stage,
or all conditional matrices if no condition is specified.

# Arguments
- `lag_matrix`: The expanded historical matrix of size (n²×n²)
- `n_stages`: Number of ahistorical stages
- `condition_stage`: If specified, return only the conditional matrix for this
  previous stage. If `nothing`, return all n conditional matrices.

# Returns
- Single condition: n×n matrix conditioned on that previous stage
- All conditions: Vector of n matrices, one per previous stage
"""
function conditional_ahistorical(lag_matrix::AbstractMatrix, n_stages::Int;
                                 condition_stage::Union{Int,Nothing}=nothing)
    n = n_stages
    expected_size = n * n
    size(lag_matrix, 1) == expected_size || throw(DimensionMismatch(
        "Expected $(expected_size)×$(expected_size) historical matrix for $n stages"))

    if condition_stage !== nothing
        return _extract_conditional(lag_matrix, n, condition_stage)
    else
        return [_extract_conditional(lag_matrix, n, s) for s in 1:n]
    end
end

"""Extract single conditional matrix for a given previous stage."""
function _extract_conditional(H::AbstractMatrix, n::Int, prev_stage::Int)
    # In the historical matrix, rows are indexed as (current_stage, next_stage)
    # and columns as (prev_stage, current_stage)
    # The conditional matrix for prev_stage s gives transitions
    # from current_stage → next_stage given previous was s
    A = zeros(eltype(H), n, n)
    for curr in 1:n
        col_idx = (prev_stage - 1) * n + curr  # column for (prev_stage, curr)
        for next in 1:n
            # Sum over all valid row entries for this next stage
            row_idx = (curr - 1) * n + next    # row for (curr, next)
            A[next, curr] = H[row_idx, col_idx]
        end
    end
    return A
end

"""
    conditional_difference(lag_matrix::AbstractMatrix, n_stages::Int)

Compute pairwise differences between conditional ahistorical matrices.

Returns an n×n matrix where entry (i,j) is the Frobenius norm of the
difference between the conditional matrix given previous stage i and
given previous stage j. High values indicate strong history-dependence
between those conditioning states.

# Arguments
- `lag_matrix`: Historical (expanded) projection matrix
- `n_stages`: Number of ahistorical stages
"""
function conditional_difference(lag_matrix::AbstractMatrix, n_stages::Int)
    cond_mats = conditional_ahistorical(lag_matrix, n_stages)
    n = n_stages
    diffs = zeros(Float64, n, n)
    for i in 1:n, j in 1:n
        diffs[i, j] = norm(cond_mats[i] .- cond_mats[j])
    end
    return diffs
end

"""
    history_dependence(lag_matrix::AbstractMatrix, n_stages::Int)

Quantify overall history-dependence of a historical MPM.

Returns the mean Frobenius distance between conditional matrices,
normalized by the mean matrix norm. Value of 0 means no history dependence;
larger values indicate stronger memory effects.
"""
function history_dependence(lag_matrix::AbstractMatrix, n_stages::Int)
    diffs = conditional_difference(lag_matrix, n_stages)
    n = n_stages
    # Mean off-diagonal difference
    total = 0.0
    count = 0
    for i in 1:n, j in 1:n
        i == j && continue
        total += diffs[i, j]
        count += 1
    end
    count == 0 && return 0.0

    # Normalize by mean matrix norm
    cond_mats = conditional_ahistorical(lag_matrix, n_stages)
    mean_norm = mean(norm.(cond_mats))
    mean_norm < eps(Float64) && return 0.0

    return (total / count) / mean_norm
end
