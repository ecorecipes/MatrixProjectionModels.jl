"""
MPMProblem — Problem type for matrix population model simulations.

Mirrors IPMProblem from IntegralProjectionModels.jl.
"""

struct MPMProblem{S<:AbstractMPMStructure, D<:AbstractDensityDependence,
                  T<:AbstractStochasticity, M, U, P, E}
    structure::S
    density::D
    stochasticity::T
    matrix::M               # MatrixProjectionModel, Vector{MPM}, or function
    n0::U                   # Initial population state
    tspan::Tuple{Int,Int}   # (t0, tf)
    p::P                    # Parameters (for density-dependent or param-resampled)
    env_state::E            # Environmental state function (for param-resampled)
    normalize::Bool          # Normalize population at each step
end

# Master constructor
function MPMProblem(structure::AbstractMPMStructure,
                    density::AbstractDensityDependence,
                    stochasticity::AbstractStochasticity,
                    matrix, n0, tspan;
                    p=nothing, env_state=nothing, normalize=false)
    MPMProblem(structure, density, stochasticity,
               matrix, n0, tspan, p, env_state, normalize)
end

# Convenience: bare matrix → LeslieMPM + DI + Det
function MPMProblem(matrix::AbstractMatrix, n0::AbstractVector, tspan::Tuple{Int,Int}; kwargs...)
    MPMProblem(LeslieMPM(), DensityIndependent(), Deterministic(),
               matrix, n0, tspan; kwargs...)
end

# Convenience: MPM → detect structure
function MPMProblem(mpm::MatrixProjectionModel, n0::AbstractVector, tspan::Tuple{Int,Int}; kwargs...)
    structure = is_leslie(mpm.A) ? LeslieMPM() : LefkovitchMPM()
    MPMProblem(structure, DensityIndependent(), Deterministic(),
               mpm, n0, tspan; kwargs...)
end

# Convenience: vector of matrices → stochastic kernel resampled
function MPMProblem(matrices::AbstractVector{<:MatrixProjectionModel}, n0::AbstractVector,
                    tspan::Tuple{Int,Int}; kwargs...)
    structure = is_leslie(matrices[1].A) ? LeslieMPM() : LefkovitchMPM()
    MPMProblem(structure, DensityIndependent(), StochasticKernelResampled(),
               matrices, n0, tspan; kwargs...)
end

# Convenience: explicit stochasticity
function MPMProblem(stoch::AbstractStochasticity, matrix, n0::AbstractVector,
                    tspan::Tuple{Int,Int}; kwargs...)
    MPMProblem(LeslieMPM(), DensityIndependent(), stoch,
               matrix, n0, tspan; kwargs...)
end

# Convenience: explicit density dependence
function MPMProblem(density::DensityDependent, matrix, n0::AbstractVector,
                    tspan::Tuple{Int,Int}; kwargs...)
    MPMProblem(LeslieMPM(), density, Deterministic(),
               matrix, n0, tspan; kwargs...)
end

# remake
function remake(prob::MPMProblem;
                structure=prob.structure, density=prob.density,
                stochasticity=prob.stochasticity, matrix=prob.matrix,
                n0=prob.n0, tspan=prob.tspan, p=prob.p,
                env_state=prob.env_state, normalize=prob.normalize)
    MPMProblem(structure, density, stochasticity,
               matrix, n0, tspan, p, env_state, normalize)
end

function Base.show(io::IO, prob::MPMProblem)
    print(io, "MPMProblem(", typeof(prob.structure).name.name, ", ",
          typeof(prob.density).name.name, ", ",
          typeof(prob.stochasticity).name.name, ", tspan=", prob.tspan, ")")
end
