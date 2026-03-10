"""
Matrix-level perturbation analysis.
Based on Rage::perturb_matrix.
"""

"""
    perturb_matrix(matA::AbstractMatrix; pert::Real=1e-6,
                   type::Symbol=:sensitivity, demog_stat::Function=lambda)

Numerical perturbation analysis of a projection matrix.
Perturbs each non-zero element by `pert` and measures change in `demog_stat`.

# Arguments
- `matA`: Projection matrix
- `pert`: Perturbation magnitude
- `type`: `:sensitivity` (absolute) or `:elasticity` (proportional)
- `demog_stat`: Demographic statistic function (default: `lambda`)

# Returns
Matrix of sensitivities or elasticities.
"""
function perturb_matrix(matA::AbstractMatrix;
                        pert::Real=1e-6,
                        type::Symbol=:sensitivity,
                        demog_stat::Function=lambda)
    n = size(matA, 1)
    base_val = demog_stat(matA)
    result = zeros(n, n)

    for j in 1:n, i in 1:n
        A_pert = copy(float.(matA))
        A_pert[i, j] += pert
        pert_val = demog_stat(A_pert)

        if type == :sensitivity
            result[i, j] = (pert_val - base_val) / pert
        elseif type == :elasticity
            if base_val != 0 && matA[i, j] != 0
                result[i, j] = ((pert_val - base_val) / pert) * (matA[i, j] / base_val)
            end
        else
            throw(ArgumentError("type must be :sensitivity or :elasticity"))
        end
    end

    return result
end
