"""
Utility functions specific to MatrixProjectionModels.jl
area_under_curve comes from StructuredPopulationCore.
"""

"""
    _validate_matrix(A::AbstractMatrix; name="A")

Validate that A is a square non-negative matrix suitable for population projection.
"""
function _validate_matrix(A::AbstractMatrix; name="A")
    size(A, 1) == size(A, 2) || throw(DimensionMismatch("$name must be square"))
    all(A .>= 0) || @warn "$name contains negative entries"
    return nothing
end

"""
    _validate_start(start::Int, n::Int)

Validate that start index is within matrix dimensions.
"""
function _validate_start(start::Int, n::Int)
    1 <= start <= n || throw(ArgumentError("start=$start must be between 1 and $n"))
end

"""
    _fundamental_matrix(U::AbstractMatrix)

Compute the fundamental matrix N = (I - U)^{-1}.
N[i,j] gives expected time spent in stage j starting from stage i.
Returns nothing if (I - U) is singular.
"""
function _fundamental_matrix(U::AbstractMatrix)
    n = size(U, 1)
    ImU = Matrix{Float64}(I, n, n) .- U
    d = det(ImU)
    abs(d) < 1e-15 && return nothing
    return inv(ImU)
end
