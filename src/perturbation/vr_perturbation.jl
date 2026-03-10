"""
Vital rate perturbation analysis.
Based on Rage::perturb_vr.
"""

"""
    perturb_vr(matU::AbstractMatrix, matF::AbstractMatrix;
               pert::Real=1e-6, type::Symbol=:sensitivity,
               demog_stat::Function=lambda,
               matC::Union{Nothing,AbstractMatrix}=nothing)

Perturbation analysis at the vital rate level.
Perturbs survival, growth, shrinkage, stasis, fecundity, and clonality
independently and measures change in demographic statistic.

# Returns
NamedTuple with fields: survival, growth, shrinkage, stasis, fecundity, clonality.
Each is a scalar (mean sensitivity/elasticity across stages).
"""
function perturb_vr(matU::AbstractMatrix, matF::AbstractMatrix;
                    pert::Real=1e-6,
                    type::Symbol=:sensitivity,
                    demog_stat::Function=lambda,
                    matC::Union{Nothing,AbstractMatrix}=nothing)
    n = size(matU, 1)
    C = matC === nothing ? zeros(n, n) : matC
    A = matU .+ matF .+ C
    base_val = demog_stat(A)

    function _perturb_component(mat_component, mat_base)
        total = 0.0
        count = 0
        for j in 1:n, i in 1:n
            if mat_component[i, j] > 0
                A_pert = copy(float.(mat_base))
                A_pert[i, j] += pert
                pert_val = demog_stat(A_pert)
                sens = (pert_val - base_val) / pert
                if type == :elasticity && base_val != 0 && mat_component[i, j] != 0
                    sens = sens * mat_component[i, j] / base_val
                end
                total += sens
                count += 1
            end
        end
        return count > 0 ? total / count : 0.0
    end

    # Decompose U into survival components
    # Stasis: diagonal of U
    stasis_mat = diagm(diag(matU))
    # Growth: upper triangle of U (transitions to earlier/smaller stages)
    growth_mat = zeros(n, n)
    for j in 1:n, i in 1:(j-1)
        growth_mat[i, j] = matU[i, j]
    end
    # Shrinkage: lower triangle of U (below diagonal)
    shrinkage_mat = zeros(n, n)
    for j in 1:n, i in (j+1):n
        shrinkage_mat[i, j] = matU[i, j]
    end

    return (
        survival = _perturb_component(matU, A),
        growth = _perturb_component(growth_mat, A),
        shrinkage = _perturb_component(shrinkage_mat, A),
        stasis = _perturb_component(stasis_mat, A),
        fecundity = _perturb_component(matF, A),
        clonality = _perturb_component(C, A)
    )
end
