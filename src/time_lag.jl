"""
Time-lagged matrix projection models.

Wraps a `MatrixProjectionModel` with lag structure information and provides
the augmented matrix for eigenanalysis and iteration.
"""

"""
    LaggedMPM{T<:Real}

A matrix projection model with time-lagged transitions.

The standard model `n(t+1) = A·n(t)` is generalized to
`n(t+1) = U·n(t) + F·n(t-F_lag) + C·n(t-C_lag)`, solved via state augmentation.

# Fields
- `mpm::MatrixProjectionModel{T, Matrix{T}}`: Original (unlagged) model
- `augmented::Matrix{T}`: `(L+1)n × (L+1)n` augmented block matrix
- `lag_structure::TimeLagStructure`: Lag specification
- `U_lag::Int`: Lag for survival/growth (default 0 = immediate)
- `F_lag::Int`: Lag for sexual fecundity (default 1)
- `C_lag::Int`: Lag for clonal fecundity (default 1)
"""
struct LaggedMPM{T<:Real}
    mpm::MatrixProjectionModel{T, Matrix{T}}
    augmented::Matrix{T}
    lag_structure::TimeLagStructure
    U_lag::Int
    F_lag::Int
    C_lag::Int
end

"""
    LaggedMPM(mpm::MatrixProjectionModel; U_lag=0, F_lag=1, C_lag=1)

Create a lagged MPM. By default, survival/growth (U) acts on the current state
and fecundity (F, C) acts on the previous time step's state.

The maximum lag determines the augmented state dimension: `(max_lag + 1) × n_stages`.
"""
function LaggedMPM(mpm::MatrixProjectionModel{T}; U_lag::Int=0, F_lag::Int=1, C_lag::Int=1) where {T}
    max_lag = max(U_lag, F_lag, C_lag)
    max_lag > 0 || throw(ArgumentError("At least one lag must be > 0"))
    lag_structure = TimeLagStructure(max_lag)

    n = n_stages(mpm)
    n_total = (max_lag + 1) * n

    # Build lag kernel vector: K_k for each lag index 0, 1, ..., max_lag
    lag_kernels = [zeros(T, n, n) for _ in 0:max_lag]
    lag_kernels[U_lag + 1] .+= mpm.U
    lag_kernels[F_lag + 1] .+= mpm.F
    lag_kernels[C_lag + 1] .+= mpm.C

    augmented = Matrix{T}(expand_lag_matrix(lag_kernels, lag_structure))

    return LaggedMPM{T}(mpm, augmented, lag_structure, U_lag, F_lag, C_lag)
end

"""
    LaggedMPM(lag_kernels::AbstractVector{<:AbstractMatrix}; lag_structure=nothing)

Create a lagged MPM directly from a vector of lag kernel matrices `[K_0, K_1, ..., K_L]`.
"""
function LaggedMPM(lag_kernels::AbstractVector{<:AbstractMatrix{T}};
        lag_structure::Union{Nothing, TimeLagStructure}=nothing) where {T}
    L = length(lag_kernels) - 1
    if lag_structure === nothing
        lag_structure = TimeLagStructure(L)
    end
    A_sum = sum(lag_kernels)
    n = size(lag_kernels[1], 1)
    mpm = MatrixProjectionModel(Matrix{T}(A_sum))
    augmented = Matrix{T}(expand_lag_matrix(lag_kernels, lag_structure))
    return LaggedMPM{T}(mpm, augmented, lag_structure, 0, min(1, L), min(1, L))
end

Base.size(lm::LaggedMPM) = size(lm.augmented)
Base.eltype(::LaggedMPM{T}) where {T} = T

function Base.show(io::IO, lm::LaggedMPM)
    n = n_stages(lm.mpm)
    L = lm.lag_structure.max_lag
    print(io, "LaggedMPM($(n) stages, max_lag=$L)")
end

# Analysis dispatches — extend ProjectionModels functions
ProjectionModels.lambda(lm::LaggedMPM) = lambda(lm.augmented)
ProjectionModels.stable_distribution(lm::LaggedMPM) = stable_distribution(lm.augmented)
ProjectionModels.reproductive_value(lm::LaggedMPM) = reproductive_value(lm.augmented)
ProjectionModels.sensitivity(lm::LaggedMPM) = sensitivity(lm.augmented)
ProjectionModels.elasticity(lm::LaggedMPM) = elasticity(lm.augmented)
ProjectionModels.damping_ratio(lm::LaggedMPM) = damping_ratio(lm.augmented)

"""
    net_repro_rate(lm::LaggedMPM; kwargs...)

Net reproductive rate for a lagged MPM, using the augmented fundamental matrix approach.
"""
function net_repro_rate(lm::LaggedMPM; kwargs...)
    result = extract_lag_components(lm.augmented, n_stages(lm.mpm), lm.lag_structure)
    net_repro_rate_lagged(result.kernels, lm.lag_structure)
end

# --- MPMProblem convenience constructor for LaggedMPM ---

function MPMProblem(lm::LaggedMPM, n0::AbstractVector, tspan::Tuple{Int,Int}; kwargs...)
    structure = is_leslie(lm.mpm.A) ? LeslieMPM() : LefkovitchMPM()
    MPMProblem(structure, DensityIndependent(), Deterministic(),
               lm, n0, tspan; kwargs...)
end

# --- Solve dispatches for LaggedMPM ---

function _get_matrix(lm::LaggedMPM)
    return lm.augmented
end

_is_lagged_matrix(::LaggedMPM) = true

function _solve_lagged(prob::MPMProblem)
    lm = prob.matrix::LaggedMPM
    n = n_stages(lm.mpm)
    L = lm.lag_structure.max_lag
    t0, tf = prob.tspan
    n_steps = tf - t0

    # Extract lag component matrices
    components = extract_lag_components(lm.augmented, n, lm.lag_structure)
    lag_kernels = components.kernels  # [K_0, K_1, ..., K_L]

    # Initialize history: all lag slots start with n0
    history = [copy(float.(prob.n0)) for _ in 0:L]

    # Store physical state (first n elements) at each time
    u = Vector{Vector{Float64}}(undef, n_steps + 1)
    u[1] = copy(float.(prob.n0))
    lambdas = Float64[]

    for t in 1:n_steps
        # n_new = Σ_k K_k · history[k]
        n_new = lag_kernels[1] * history[1]
        for k in 1:L
            n_new .+= lag_kernels[k + 1] * history[k + 1]
        end

        pop_t = sum(history[1])
        pop_new = sum(n_new)
        push!(lambdas, pop_t > 0 ? pop_new / pop_t : 0.0)

        if prob.normalize && pop_new > 0
            n_new ./= pop_new
        end

        # Shift history: drop oldest, insert new
        for k in L:-1:1
            history[k + 1] = history[k]
        end
        history[1] = n_new

        u[t + 1] = copy(n_new)
    end

    ts = collect(t0:tf)
    ea = eigenanalysis_power(lm.augmented)
    return MPMSolution(ts, u, lm.augmented, ea, :Success, lambdas)
end
