"""
MPM splitting: decompose A into U + F + C.
Based on Rage::mpm_split.
"""

"""
    mpm_split(matA::AbstractMatrix; matC::Union{Nothing,AbstractMatrix}=nothing)

Attempt to decompose a projection matrix A into U (survival) + F (fecundity).
Uses a simple heuristic: first row = F, rest = U.

If matC is provided, A = U + F + C.

# Returns
`MatrixProjectionModel` with decomposition.
"""
function mpm_split(matA::AbstractMatrix; matC::Union{Nothing,AbstractMatrix}=nothing)
    n = size(matA, 1)
    U = copy(float.(matA))
    F = zeros(Float64, n, n)
    C = matC === nothing ? zeros(Float64, n, n) : copy(float.(matC))

    # Heuristic: first row is fecundity
    F[1, :] .= U[1, :]
    U[1, :] .= 0.0

    # Subtract clonal if provided
    if matC !== nothing
        U .-= C
    end

    return MatrixProjectionModel(U, F, C)
end
