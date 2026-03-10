"""
Matrix properties specific to MatrixProjectionModels.
is_irreducible, is_primitive, is_ergodic come from ProjectionModels.
"""

"""
    is_leslie(A::AbstractMatrix)

Test if matrix A has Leslie structure:
- Non-negative entries in first row (fecundity)
- Positive entries on sub-diagonal (survival)
- All other entries are zero
"""
function is_leslie(A::AbstractMatrix)
    n = size(A, 1)
    n < 2 && return false
    for i in 1:n, j in 1:n
        if i == 1
            # First row: must be non-negative
            A[i, j] < 0 && return false
        elseif i == j + 1
            # Sub-diagonal: must be positive
            A[i, j] <= 0 && return false
        else
            # All other entries: must be zero
            A[i, j] != 0 && return false
        end
    end
    return true
end

# MatrixProjectionModel dispatches — extend ProjectionModels functions
ProjectionModels.is_irreducible(m::MatrixProjectionModel) = is_irreducible(m.A)
ProjectionModels.is_primitive(m::MatrixProjectionModel) = is_primitive(m.A)
ProjectionModels.is_ergodic(m::MatrixProjectionModel) = is_ergodic(m.A)
is_leslie(m::MatrixProjectionModel) = is_leslie(m.A)
