"""
MPM collapse: reduce stages by merging.
Based on Rage::mpm_collapse.
"""

"""
    mpm_collapse(matU::AbstractMatrix, matF::AbstractMatrix,
                 collapse::Vector{Vector{Int}};
                 matC::Union{Nothing,AbstractMatrix}=nothing)

Collapse stages of an MPM according to a grouping specification.

# Arguments
- `matU`: Survival matrix
- `matF`: Fecundity matrix
- `collapse`: Vector of vectors, each containing indices of stages to merge.
  E.g., `[[1,2], [3], [4,5]]` merges stages 1&2 and 4&5.
- `matC`: Optional clonal matrix

# Returns
`MatrixProjectionModel` with collapsed stages.

# Method
Uses the stable stage distribution as weights for collapsing.
The dominant right eigenvector of A determines the relative contribution
of each stage within a group.
"""
function mpm_collapse(matU::AbstractMatrix, matF::AbstractMatrix,
                      collapse::Vector{Vector{Int}};
                      matC::Union{Nothing,AbstractMatrix}=nothing)
    n = size(matU, 1)
    C = matC === nothing ? zeros(n, n) : matC
    A = matU .+ matF .+ C
    n_new = length(collapse)

    # Get stable stage distribution for weighting
    w = stable_distribution(A)
    w = max.(w, 0.0)

    # Build collapse mapping matrix T (n_new × n)
    # T[k, i] = weight of old stage i in new stage k
    T_mat = zeros(n_new, n)
    for (k, group) in enumerate(collapse)
        group_total = sum(w[i] for i in group)
        for i in group
            T_mat[k, i] = group_total > 0 ? w[i] / group_total : 1.0 / length(group)
        end
    end

    # Build expansion matrix S (n × n_new)
    # S[i, k] = 1 if old stage i belongs to new stage k
    S_mat = zeros(n, n_new)
    for (k, group) in enumerate(collapse)
        for i in group
            S_mat[i, k] = 1.0
        end
    end

    # Collapsed matrices: M_new = T * M * S
    U_new = T_mat * matU * S_mat
    F_new = T_mat * matF * S_mat
    C_new = T_mat * C * S_mat

    stage_names = [Symbol("stage_$i") for i in 1:n_new]
    return MatrixProjectionModel(U_new, F_new, C_new; stage_names=stage_names)
end
