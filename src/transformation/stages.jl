"""
Stage classification utilities.
Based on Rage::repro_stages, standard_stages.
"""

"""
    repro_stages(matR::AbstractMatrix)

Identify reproductive stages: stages that produce offspring.
A stage j is reproductive if any element in column j of matR is positive.

# Returns
`BitVector` of length n.
"""
function repro_stages(matR::AbstractMatrix)
    n = size(matR, 2)
    return BitVector([any(matR[:, j] .> 0) for j in 1:n])
end

"""
    standard_stages(matF::AbstractMatrix; matC::Union{Nothing,AbstractMatrix}=nothing)

Classify stages into standard types based on matrix structure.

# Returns
Vector of `StageClassType`: ActiveStage, PropaguleStage, or DormantStage.
"""
function standard_stages(matF::AbstractMatrix;
                         matC::Union{Nothing,AbstractMatrix}=nothing)
    n = size(matF, 1)
    C = matC === nothing ? zeros(n, n) : matC

    stages = fill(ActiveStage, n)

    # Propagule stages: stages that receive reproduction but don't themselves reproduce
    repro = repro_stages(matF .+ C)
    receives_repro = BitVector([any((matF .+ C)[i, :] .> 0) for i in 1:n])

    for i in 1:n
        if receives_repro[i] && !repro[i]
            stages[i] = PropaguleStage
        end
    end

    return stages
end

"""
    name_stages(mpm::MatrixProjectionModel, names::AbstractVector{<:AbstractString})

Create a copy of the MPM with updated stage names.
"""
function name_stages(mpm::MatrixProjectionModel, names::AbstractVector{<:AbstractString})
    length(names) == n_stages(mpm) ||
        throw(ArgumentError("names must have length $(n_stages(mpm))"))
    stage_names = Symbol.(names)
    return MatrixProjectionModel(mpm.A, mpm.U, mpm.F, mpm.C, mpm.stages, stage_names)
end
