"""
MPM standardization and rearrangement.
Based on Rage::mpm_standardize, mpm_rearrange.
"""

"""
    mpm_standardize(matU::AbstractMatrix, matF::AbstractMatrix;
                    repro_stages::Union{Nothing,AbstractVector{Bool}}=nothing,
                    matC::Union{Nothing,AbstractMatrix}=nothing)

Standardize an MPM so that:
1. The first stage is a propagule/newborn stage
2. Reproductive stages produce into the first stage only

If repro_stages is not provided, it's detected from matF.

# Returns
`MatrixProjectionModel` with standardized structure.
"""
function mpm_standardize(matU::AbstractMatrix, matF::AbstractMatrix;
                         repro_stages::Union{Nothing,AbstractVector{Bool}}=nothing,
                         matC::Union{Nothing,AbstractMatrix}=nothing)
    n = size(matU, 1)
    C = matC === nothing ? zeros(n, n) : copy(float.(matC))

    if repro_stages === nothing
        repro_stages = MatrixProjectionModels.repro_stages(matF)
    end

    U_new = copy(float.(matU))
    F_new = zeros(Float64, n, n)

    # Standardize: all reproduction goes to first row
    for j in 1:n
        total_repro = sum(matF[:, j])
        F_new[1, j] = total_repro
    end

    return MatrixProjectionModel(U_new, F_new, C)
end

"""
    mpm_rearrange(matU::AbstractMatrix, matF::AbstractMatrix;
                  new_order::AbstractVector{Int},
                  matC::Union{Nothing,AbstractMatrix}=nothing)

Rearrange stages of an MPM according to a new ordering.

# Arguments
- `new_order`: Permutation vector. `new_order[k] = i` means old stage i becomes new stage k.
"""
function mpm_rearrange(matU::AbstractMatrix, matF::AbstractMatrix;
                       new_order::AbstractVector{Int},
                       matC::Union{Nothing,AbstractMatrix}=nothing)
    n = size(matU, 1)
    length(new_order) == n || throw(ArgumentError("new_order must have length $n"))
    sort(new_order) == 1:n || throw(ArgumentError("new_order must be a permutation of 1:$n"))

    C = matC === nothing ? zeros(n, n) : matC

    U_new = matU[new_order, new_order]
    F_new = matF[new_order, new_order]
    C_new = C[new_order, new_order]

    return MatrixProjectionModel(float.(U_new), float.(F_new), float.(C_new))
end
