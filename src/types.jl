"""
Core types for MatrixProjectionModels.jl

MatrixProjectionModel: A = U + F + C decomposition
Trait types for dispatch: structure, density dependence, stochasticity, algorithms
"""

# --- Stage Classification ---

@enum StageClassType ActiveStage PropaguleStage DormantStage

struct StageClass
    organized::StageClassType
    author::String
    number::Int
end

# --- Matrix Population Model ---

"""
    MatrixProjectionModel{T<:Real, M<:AbstractMatrix{T}}

A matrix population model with A = U + F + C decomposition.
Parameterized on both element type `T` and matrix type `M`,
so you can use `Matrix`, `SMatrix`, `SparseMatrixCSC`, etc.

# Fields
- `A::M`: Full projection matrix
- `U::M`: Survival/growth transitions
- `F::M`: Sexual reproduction
- `C::M`: Clonal reproduction
- `stages::Vector{StageClass}`: Stage classification metadata
- `stage_names::Vector{Symbol}`: Names for each stage
"""
struct MatrixProjectionModel{T<:Real, M<:AbstractMatrix{T}}
    A::M
    U::M
    F::M
    C::M
    stages::Vector{StageClass}
    stage_names::Vector{Symbol}

    function MatrixProjectionModel{T,M}(A, U, F, C, stages, stage_names) where {T, M}
        n = size(A, 1)
        size(A, 2) == n || throw(DimensionMismatch("A must be square"))
        size(U) == (n, n) || throw(DimensionMismatch("U must match A dimensions"))
        size(F) == (n, n) || throw(DimensionMismatch("F must match A dimensions"))
        size(C) == (n, n) || throw(DimensionMismatch("C must match A dimensions"))
        length(stage_names) == n || throw(DimensionMismatch("stage_names length must match matrix dimension"))
        new{T,M}(A, U, F, C, stages, stage_names)
    end
end

# Infer M from the matrix arguments
function _mpm_construct(A::MA, U::MU, F::MF, C::MC,
                        stages, stage_names) where {T, MA<:AbstractMatrix{T},
                                                    MU<:AbstractMatrix{T},
                                                    MF<:AbstractMatrix{T},
                                                    MC<:AbstractMatrix{T}}
    # Use the type of A as the canonical matrix type
    MatrixProjectionModel{T,MA}(A, U, F, C, stages, stage_names)
end

function _default_stage_names(n)
    [Symbol("stage_$i") for i in 1:n]
end

# Constructor: full specification (positional stages + stage_names)
function MatrixProjectionModel(A::MA, U::MU, F::MF, C::MC,
                               stages::Vector{StageClass},
                               stage_names::Vector{Symbol}) where {T,
                                   MA<:AbstractMatrix{T}, MU<:AbstractMatrix{T},
                                   MF<:AbstractMatrix{T}, MC<:AbstractMatrix{T}}
    _mpm_construct(A, U, F, C, stages, stage_names)
end

# Constructor: from A + U + F + C (use provided A, don't recompute)
function MatrixProjectionModel(A::MA, U::MU, F::MF, C::MC;
                               stages::Vector{StageClass}=StageClass[],
                               stage_names::Vector{Symbol}=Symbol[]) where {T,
                                   MA<:AbstractMatrix{T}, MU<:AbstractMatrix{T},
                                   MF<:AbstractMatrix{T}, MC<:AbstractMatrix{T}}
    n = size(A, 1)
    if isempty(stage_names)
        stage_names = _default_stage_names(n)
    end
    _mpm_construct(A, U, F, C, stages, stage_names)
end

# Constructor: from U + F + C (compute A = U + F + C)
function MatrixProjectionModel(U::MU, F::MF, C::MC;
                               stages::Vector{StageClass}=StageClass[],
                               stage_names::Vector{Symbol}=Symbol[]) where {T,
                                   MU<:AbstractMatrix{T}, MF<:AbstractMatrix{T},
                                   MC<:AbstractMatrix{T}}
    A = U .+ F .+ C
    MatrixProjectionModel(A, U, F, C; stages=stages, stage_names=stage_names)
end

# Constructor: from U + F (no clonal, C = 0)
function MatrixProjectionModel(U::MU, F::MF;
                               stages::Vector{StageClass}=StageClass[],
                               stage_names::Vector{Symbol}=Symbol[]) where {T,
                                   MU<:AbstractMatrix{T}, MF<:AbstractMatrix{T}}
    C = zero(U)
    MatrixProjectionModel(U, F, C; stages=stages, stage_names=stage_names)
end

# Constructor: from A only (U=A, F=0, C=0 — unknown decomposition)
function MatrixProjectionModel(A::MA;
                               stages::Vector{StageClass}=StageClass[],
                               stage_names::Vector{Symbol}=Symbol[]) where {T,
                                   MA<:AbstractMatrix{T}}
    U = copy(A)
    F = zero(A)
    C = zero(A)
    n = size(A, 1)
    if isempty(stage_names)
        stage_names = _default_stage_names(n)
    end
    _mpm_construct(A, U, F, C, stages, stage_names)
end

# Constructor: from survival and fecundity vectors (Leslie matrix → always Matrix{T})
function MatrixProjectionModel(survival::AbstractVector{T}, fecundity::AbstractVector{T};
                               stages::Vector{StageClass}=StageClass[],
                               stage_names::Vector{Symbol}=Symbol[]) where {T<:Real}
    mpm = make_leslie_mpm(survival, fecundity)
    resolved_stage_names = isempty(stage_names) ? mpm.stage_names : stage_names
    return MatrixProjectionModel(mpm.A, mpm.U, mpm.F, mpm.C;
                                 stages=stages, stage_names=resolved_stage_names)
end

# Type promotion constructors (produce Matrix{T})
function MatrixProjectionModel(A::AbstractMatrix, U::AbstractMatrix,
                               F::AbstractMatrix, C::AbstractMatrix; kwargs...)
    T = promote_type(eltype(A), eltype(U), eltype(F), eltype(C))
    MatrixProjectionModel(Matrix{T}(A), Matrix{T}(U), Matrix{T}(F), Matrix{T}(C); kwargs...)
end

function MatrixProjectionModel(U::AbstractMatrix, F::AbstractMatrix,
                               C::AbstractMatrix; kwargs...)
    T = promote_type(eltype(U), eltype(F), eltype(C))
    MatrixProjectionModel(Matrix{T}(U), Matrix{T}(F), Matrix{T}(C); kwargs...)
end

function MatrixProjectionModel(U::AbstractMatrix, F::AbstractMatrix; kwargs...)
    T = promote_type(eltype(U), eltype(F))
    MatrixProjectionModel(Matrix{T}(U), Matrix{T}(F); kwargs...)
end

# --- AbstractMatrix interface forwarding to A ---

Base.size(m::MatrixProjectionModel) = size(m.A)
Base.size(m::MatrixProjectionModel, d::Integer) = size(m.A, d)
Base.getindex(m::MatrixProjectionModel, i...) = getindex(m.A, i...)
Base.eltype(::Type{<:MatrixProjectionModel{T}}) where {T} = T
Base.eltype(m::MatrixProjectionModel) = eltype(m.A)
Base.length(m::MatrixProjectionModel) = length(m.A)
Base.iterate(m::MatrixProjectionModel, state...) = iterate(m.A, state...)

# Allow Matrix(mpm) conversion
Base.Matrix(m::MatrixProjectionModel) = Matrix(m.A)
Base.convert(::Type{Matrix{T}}, m::MatrixProjectionModel) where {T} = convert(Matrix{T}, m.A)

# Number of stages
n_stages(m::MatrixProjectionModel) = size(m.A, 1)

function Base.show(io::IO, m::MatrixProjectionModel{T,M}) where {T,M}
    n = n_stages(m)
    mname = M === Matrix{T} ? "$T" : "$T, $M"
    print(io, "MatrixProjectionModel{$mname} ($n stages)")
end

function Base.show(io::IO, ::MIME"text/plain", m::MatrixProjectionModel{T,M}) where {T,M}
    n = n_stages(m)
    mname = M === Matrix{T} ? "$T" : "$T, $M"
    println(io, "MatrixProjectionModel{$mname} with $n stages:")
    println(io, "  Stage names: ", m.stage_names)
    println(io, "  A (projection matrix):")
    show(io, "text/plain", m.A)
end

# --- Trait Types ---

# Structure (MPM-specific)
abstract type AbstractMPMStructure <: AbstractProjectionStructure end
struct LeslieMPM <: AbstractMPMStructure end
struct LefkovitchMPM <: AbstractMPMStructure end

# Density dependence, stochasticity, and algorithm types are
# imported from StructuredPopulationCore via `using StructuredPopulationCore`
