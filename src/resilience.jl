# Resilience operators: `ResilienceCore.linearise` for matrix projection models.

"""
    linearise(m::MatrixProjectionModel) -> ResilienceCore.ProjectionMatrix
    linearise(prob::MPMProblem)         -> ResilienceCore.ProjectionMatrix

The projection operator of a density-independent, deterministic matrix projection model, for
the stability and transient metrics of [ResilienceCore.jl](https://github.com/ecorecipes/ResilienceCore.jl)
(`growth_rate`, `resilience`, `reactivity`, `max_amplification`, …).

Only such models have a single projection matrix. A density-dependent `MPMProblem` is linearised
at a fixed point instead, with `linearise(to_discrete_problem(prob), n_star)`; a stochastic one has
a stochastic growth rate but no single spectrum. Both throw an `ArgumentError`.
"""
ResilienceCore.linearise(m::MatrixProjectionModel) = ResilienceCore.ProjectionMatrix(m.A)

function ResilienceCore.linearise(prob::MPMProblem)
    prob.density isa StructuredPopulationCore.DensityIndependent ||
        throw(ArgumentError("a density-dependent MPM has no single projection matrix; linearise its map at a " *
                            "fixed point instead, with `linearise(to_discrete_problem(prob), n_star)`"))
    prob.stochasticity isa StructuredPopulationCore.Deterministic ||
        throw(ArgumentError("a stochastic MPM has no single projection matrix; use its stochastic growth rate"))
    M = prob.matrix
    M isa MatrixProjectionModel && return ResilienceCore.ProjectionMatrix(M.A)
    M isa AbstractMatrix && return ResilienceCore.ProjectionMatrix(M)
    throw(ArgumentError("cannot read a projection matrix from an MPMProblem holding a $(typeof(M))"))
end
