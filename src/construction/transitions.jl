"""
Sparse transition constructors for MatrixProjectionModel.

Allows specifying projection matrices via `(from => to) => value` pairs
instead of writing full matrices.

# Examples
```julia
# A-only from sparse transitions
mpm = MatrixProjectionModel([:seed, :small, :large],
    (:seed => :small) => 0.2,
    (:small => :large) => 0.4,
    (:small => :small) => 0.3,
    (:large => :large) => 0.7,
    (:large => :seed) => 5.0)

# With U/F/C decomposition
mpm = MatrixProjectionModel([:seed, :small, :large];
    U = [(:seed => :small) => 0.2, (:small => :large) => 0.4,
         (:small => :small) => 0.3, (:large => :large) => 0.7],
    F = [(:large => :seed) => 5.0])
```
"""

"""
    _transitions_to_matrix(T, stage_idx, n, entries)

Build an n×n matrix from sparse `(from => to) => value` entries.
Convention: `(from => to) => value` sets `M[to_idx, from_idx] = value`.
Multiple entries to the same cell are summed.
"""
function _transitions_to_matrix(::Type{T}, stage_idx::Dict{Symbol,Int}, n::Int,
                                entries) where {T<:Real}
    M = zeros(T, n, n)
    for (pair, val) in entries
        from, to = pair
        haskey(stage_idx, from) || throw(ArgumentError("Unknown stage name :$from"))
        haskey(stage_idx, to) || throw(ArgumentError("Unknown stage name :$to"))
        M[stage_idx[to], stage_idx[from]] += val
    end
    return M
end

# Constructor: A-only from sparse transitions (U=A, F=0, C=0)
"""
    MatrixProjectionModel(stage_names, transitions...; stages=StageClass[])

Construct a `MatrixProjectionModel` from named stages and sparse
`(from => to) => value` transition pairs. Sets A from the transitions,
with U=A, F=0, C=0 (unknown decomposition).
"""
function MatrixProjectionModel(stage_names::AbstractVector{Symbol},
                               transitions::Pair{Pair{Symbol,Symbol}, <:Real}...;
                               stages::Vector{StageClass}=StageClass[])
    n = length(stage_names)
    stage_idx = Dict(s => i for (i, s) in enumerate(stage_names))
    T = isempty(transitions) ? Float64 : promote_type(map(t -> typeof(t.second), transitions)...)
    T = T <: AbstractFloat ? T : Float64
    A = _transitions_to_matrix(T, stage_idx, n, transitions)
    U = copy(A)
    F = zeros(T, n, n)
    C = zeros(T, n, n)
    _mpm_construct(A, U, F, C, stages, collect(Symbol, stage_names))
end

# Constructor: U/F/C from sparse transitions
"""
    MatrixProjectionModel(stage_names; U=[], F=[], C=[], stages=StageClass[])

Construct a `MatrixProjectionModel` from named stages and sparse
`(from => to) => value` entries for U, F, and C matrices.
Computes A = U + F + C.
"""
function MatrixProjectionModel(stage_names::AbstractVector{Symbol};
                               U::AbstractVector=Pair{Pair{Symbol,Symbol},Float64}[],
                               F::AbstractVector=Pair{Pair{Symbol,Symbol},Float64}[],
                               C::AbstractVector=Pair{Pair{Symbol,Symbol},Float64}[],
                               stages::Vector{StageClass}=StageClass[])
    n = length(stage_names)
    stage_idx = Dict(s => i for (i, s) in enumerate(stage_names))

    # Determine element type from all entries
    all_vals = vcat(
        [t.second for t in U],
        [t.second for t in F],
        [t.second for t in C]
    )
    T = isempty(all_vals) ? Float64 : promote_type(map(typeof, all_vals)...)
    T = T <: AbstractFloat ? T : Float64

    U_mat = _transitions_to_matrix(T, stage_idx, n, U)
    F_mat = _transitions_to_matrix(T, stage_idx, n, F)
    C_mat = _transitions_to_matrix(T, stage_idx, n, C)
    A_mat = U_mat .+ F_mat .+ C_mat

    _mpm_construct(A_mat, U_mat, F_mat, C_mat, stages, collect(Symbol, stage_names))
end
