"""
Transition-level perturbation analysis.
Based on Rage::perturb_trans.
"""

"""
    perturb_trans(matU::AbstractMatrix, matF::AbstractMatrix;
                  pert::Real=1e-6, type::Symbol=:sensitivity,
                  demog_stat::Function=lambda,
                  matC::Union{Nothing,AbstractMatrix}=nothing)

Perturbation analysis at the transition level.
Returns separate sensitivity/elasticity matrices for U, F, and C components.

# Returns
NamedTuple (U, F, C) of perturbation matrices.
"""
function perturb_trans(matU::AbstractMatrix, matF::AbstractMatrix;
                       pert::Real=1e-6,
                       type::Symbol=:sensitivity,
                       demog_stat::Function=lambda,
                       matC::Union{Nothing,AbstractMatrix}=nothing)
    n = size(matU, 1)
    C = matC === nothing ? zeros(n, n) : matC
    A = matU .+ matF .+ C
    base_val = demog_stat(A)

    function _perturb_mat(component)
        result = zeros(n, n)
        for j in 1:n, i in 1:n
            if component[i, j] > 0 || type == :sensitivity
                A_pert = copy(float.(A))
                A_pert[i, j] += pert
                pert_val = demog_stat(A_pert)
                sens = (pert_val - base_val) / pert
                if type == :elasticity
                    if base_val != 0 && component[i, j] != 0
                        result[i, j] = sens * component[i, j] / base_val
                    end
                else
                    result[i, j] = sens
                end
            end
        end
        return result
    end

    return (U=_perturb_mat(matU), F=_perturb_mat(matF), C=_perturb_mat(C))
end
