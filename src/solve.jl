"""
CommonSolve.solve dispatches for MPMProblem → MPMSolution.
"""

struct MPMSolution{T,U,K,E} <: AbstractProjectionSolution
    t::T                         # Time steps
    u::U                         # Population states at each time step
    kernel_matrices::K           # Matrix or vector of matrices used
    eigenanalysis::E             # NamedTuple or nothing
    retcode::Symbol
    lambdas::Vector{Float64}     # Per-timestep growth rates
end

function Base.show(io::IO, sol::MPMSolution)
    nt = length(sol.t)
    ret = sol.retcode
    print(io, "MPMSolution($nt timesteps, retcode=$ret)")
end

# --- Extract matrix from problem ---

function _get_matrix(prob::MPMProblem)
    m = prob.matrix
    if m isa MatrixProjectionModel
        return m.A
    elseif m isa AbstractMatrix
        return m
    elseif hasmethod(_get_matrix, Tuple{typeof(m)})
        return _get_matrix(m)
    else
        error("Cannot extract matrix from $(typeof(m))")
    end
end

function _get_matrix(mpm::MatrixProjectionModel)
    return mpm.A
end


# --- Main solve dispatch ---

function CommonSolve.solve(prob::MPMProblem, alg::EigenAnalysis=EigenAnalysis(); kwargs...)
    A = _get_matrix(prob)
    ea = eigenanalysis_full(A)
    u = [copy(prob.n0)]
    K = prob.matrix isa MatrixProjectionModel ? prob.matrix.A : A
    return MPMSolution(collect(prob.tspan[1]:prob.tspan[1]), u, K, ea, :Success, Float64[])
end

function CommonSolve.solve(prob::MPMProblem, alg::DirectIteration; kwargs...)
    if _is_lagged_matrix(prob.matrix)
        return _solve_lagged(prob)
    end
    _solve(prob.structure, prob.density, prob.stochasticity, prob, alg; kwargs...)
end

# Fallback: non-lagged matrices
_is_lagged_matrix(m) = false
# _is_lagged_matrix(::LaggedMPM) = true is defined in time_lag.jl

# --- Deterministic, density-independent ---

function _solve(::AbstractMPMStructure, ::DensityIndependent, ::Deterministic,
                prob::MPMProblem, ::DirectIteration; kwargs...)
    A = _get_matrix(prob)
    t0, tf = prob.tspan
    n_steps = tf - t0

    u = Vector{Vector{Float64}}(undef, n_steps + 1)
    u[1] = copy(float.(prob.n0))
    lambdas = Float64[]
    u_new = similar(u[1])

    for t in 1:n_steps
        mul!(u_new, A, u[t])
        n_total = sum(u_new)
        push!(lambdas, n_total / sum(u[t]))
        if prob.normalize && n_total > 0
            u_new ./= n_total
        end
        u[t+1] = copy(u_new)
    end

    ts = collect(t0:tf)
    ea = eigenanalysis_power(A)
    return MPMSolution(ts, u, A, ea, :Success, lambdas)
end

# --- Density-dependent ---

function _solve(::AbstractMPMStructure, ::DensityDependent, ::Deterministic,
                prob::MPMProblem, ::DirectIteration; kwargs...)
    t0, tf = prob.tspan
    n_steps = tf - t0

    u = Vector{Vector{Float64}}(undef, n_steps + 1)
    u[1] = copy(float.(prob.n0))
    lambdas = Float64[]
    matrices = Matrix{Float64}[]

    for t in 1:n_steps
        # matrix is a function of (n, p, t)
        A = prob.matrix(u[t], prob.p, t0 + t - 1)
        push!(matrices, A)
        u_new = A * u[t]
        n_total = sum(u_new)
        push!(lambdas, n_total / sum(u[t]))
        if prob.normalize && n_total > 0
            u_new ./= n_total
        end
        u[t+1] = u_new
    end

    ts = collect(t0:tf)
    return MPMSolution(ts, u, matrices, nothing, :Success, lambdas)
end

# --- Stochastic kernel resampled ---

function _solve(::AbstractMPMStructure, ::DensityIndependent, ::StochasticKernelResampled,
                prob::MPMProblem, ::DirectIteration;
                rng::AbstractRNG=Random.default_rng(), kernel_seq=nothing, kwargs...)
    t0, tf = prob.tspan
    n_steps = tf - t0

    matrices = prob.matrix
    if matrices isa AbstractVector{<:MatrixProjectionModel}
        As = [m.A for m in matrices]
    else
        As = matrices
    end

    u = Vector{Vector{Float64}}(undef, n_steps + 1)
    u[1] = copy(float.(prob.n0))
    lambdas = Float64[]
    used_matrices = Matrix{Float64}[]

    for t in 1:n_steps
        idx = kernel_seq !== nothing ? kernel_seq[t] : rand(rng, 1:length(As))
        A = As[idx]
        push!(used_matrices, A)
        u_new = A * u[t]
        n_total = sum(u_new)
        push!(lambdas, n_total / sum(u[t]))
        if prob.normalize && n_total > 0
            u_new ./= n_total
        end
        u[t+1] = u_new
    end

    ts = collect(t0:tf)
    return MPMSolution(ts, u, used_matrices, nothing, :Success, lambdas)
end

# --- Stochastic parameter resampled ---

function _solve(::AbstractMPMStructure, ::DensityIndependent, ::StochasticParameterResampled,
                prob::MPMProblem, ::DirectIteration;
                rng::AbstractRNG=Random.default_rng(), kwargs...)
    t0, tf = prob.tspan
    n_steps = tf - t0

    u = Vector{Vector{Float64}}(undef, n_steps + 1)
    u[1] = copy(float.(prob.n0))
    lambdas = Float64[]
    matrices = Matrix{Float64}[]

    for t in 1:n_steps
        # env_state returns parameters for this time step
        params = prob.env_state(rng, t0 + t - 1)
        # matrix is a function of parameters
        A = prob.matrix(params)
        if A isa MatrixProjectionModel
            A = A.A
        end
        push!(matrices, A)
        u_new = A * u[t]
        n_total = sum(u_new)
        push!(lambdas, n_total / sum(u[t]))
        if prob.normalize && n_total > 0
            u_new ./= n_total
        end
        u[t+1] = u_new
    end

    ts = collect(t0:tf)
    return MPMSolution(ts, u, matrices, nothing, :Success, lambdas)
end
