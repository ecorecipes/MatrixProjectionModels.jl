"""
Stateful discrete matrix population systems.

This layer complements `MPMProblem` for models that need multiple named
components, auxiliary mutable state, scheduled events, and post-step
transition rules without hand-written outer iteration loops.

The components, the system, the built-in events and rules and the helpers that read and
change a system are defined in StructuredPopulationCore (`src/coupled_systems.jl`) and
shared with PhysiologicallyBasedDemographicModels.jl. The engine that steps a system (the
substeps, the stage fluxes and transfers that rules return, `apply_rule`, `CustomEvent` and
the loop over time steps) is also defined there (`src/coupled_stepping.jl`) and shared with
IntegralProjectionModels.jl. This file adds what is specific to matrix models: the
`MatrixProjectionModel` methods of the component hooks, the allocation pool, how a component
is turned into its projection matrix, and the problem, solution and solver.
"""

# ---------------------------------------------------------------------------
# Population components and systems
# ---------------------------------------------------------------------------

function component_stage_names(model::MatrixProjectionModel, n::Int)
    size(model.A) == (n, n) || throw(DimensionMismatch(
        "model has size $(size(model.A)); expected ($n, $n)"))
    return copy(model.stage_names)
end

function validate_component_model(model::MatrixProjectionModel, n::Int)
    size(model.A) == (n, n) || throw(DimensionMismatch(
        "model has size $(size(model.A)); expected ($n, $n)"))
end

# ---------------------------------------------------------------------------
# Hybrid allocation helpers
# ---------------------------------------------------------------------------

"""
    PriorityAllocationPool(supply, demands, labels)

Priority-ordered resource pool for hybrid daily stepping. Supply is allocated to
demands in order, with later entries receiving whatever remains.
"""
struct PriorityAllocationPool{T<:Real}
    supply::T
    demands::Vector{T}
    labels::Vector{Symbol}

    function PriorityAllocationPool(supply::T, demands::Vector{T}, labels::Vector{Symbol}) where {T<:Real}
        length(demands) == length(labels) || throw(DimensionMismatch(
            "demands and labels must have the same length"))
        new{T}(supply, demands, labels)
    end
end

function PriorityAllocationPool(supply::Real, demands::AbstractVector{<:Real},
        labels::AbstractVector{Symbol} = Symbol[])
    d = collect(float.(demands))
    T = isempty(d) ? typeof(float(supply)) : promote_type(typeof(float(supply)), eltype(d))
    lab = isempty(labels) ? [Symbol("slot_", i) for i in eachindex(d)] : collect(labels)
    return PriorityAllocationPool(T(supply), T.(d), lab)
end

"""
    allocate(pool::PriorityAllocationPool)

Allocate supply to demands in priority order.
"""
function allocate(pool::PriorityAllocationPool{T}) where {T}
    allocations = zeros(T, length(pool.demands))
    remaining = pool.supply
    for i in eachindex(pool.demands)
        allocations[i] = min(remaining, pool.demands[i])
        remaining -= allocations[i]
        remaining <= zero(T) && break
    end
    return allocations
end

"""
    supply_demand_index(pool::PriorityAllocationPool)

Return the overall supply/demand ratio in `[0, 1]`.
"""
function supply_demand_index(pool::PriorityAllocationPool)
    total_demand = sum(pool.demands)
    total_demand <= 0 && return one(pool.supply)
    return min(one(pool.supply), pool.supply / total_demand)
end

"""
    allocation_stress(pool::PriorityAllocationPool)

Return `(allocations, stress)` where stress is the unmet-demand fraction for each
slot in `[0, 1]`.
"""
function allocation_stress(pool::PriorityAllocationPool{T}) where {T}
    allocations = allocate(pool)
    stress = similar(allocations)
    for i in eachindex(allocations)
        demand = pool.demands[i]
        stress[i] = demand > zero(T) ? max(zero(T), one(T) - allocations[i] / demand) : zero(T)
    end
    return allocations, stress
end

# ---------------------------------------------------------------------------
# Coupled/stateful MPM problem and solution
# ---------------------------------------------------------------------------

"""
    CoupledMPMProblem(system, tspan; p, rules, events, observables, normalize)

Problem type for stateful, multi-component discrete matrix population systems.
"""
struct CoupledMPMProblem{Sys, P, S, R, E, O}
    system::Sys
    tspan::Tuple{Int, Int}
    p::P
    substeps::S
    rules::R
    events::E
    observables::O
    normalize::Bool
end

const StateDependentMPMProblem = CoupledMPMProblem
const HybridMPMProblem = CoupledMPMProblem

function CoupledMPMProblem(system::PopulationSystem, tspan::Tuple{Int, Int};
        p = nothing,
        substeps::AbstractVector = AbstractDailySubstep[],
        rules::AbstractVector = AbstractTransitionRule[],
        events::AbstractVector = AbstractScheduledEvent[],
        observables::AbstractVector = Observable[],
        normalize::Bool = false)
    tspan[2] >= tspan[1] || throw(ArgumentError("tspan must satisfy tspan[2] >= tspan[1]"))
    _ensure_unique_names(substeps, "substep")
    _ensure_unique_names(rules, "rule")
    _ensure_unique_names(observables, "observable")
    _require_vector_carriers(system)
    return CoupledMPMProblem(system, tspan, p, collect(substeps), collect(rules), collect(events), collect(observables), normalize)
end

function remake(prob::CoupledMPMProblem;
        system = prob.system,
        tspan = prob.tspan,
        p = prob.p,
        substeps = prob.substeps,
        rules = prob.rules,
        events = prob.events,
        observables = prob.observables,
        normalize = prob.normalize)
    CoupledMPMProblem(system, tspan;
        p = p,
        substeps = substeps,
        rules = rules,
        events = events,
        observables = observables,
        normalize = normalize)
end

function Base.show(io::IO, prob::CoupledMPMProblem)
    print(io,
        "CoupledMPMProblem($(length(prob.system)) components, ",
        "substeps=$(length(prob.substeps)), rules=$(length(prob.rules)), ",
        "events=$(length(prob.events)), ",
        "tspan=$(prob.tspan))")
end

"""
    CoupledMPMSolution

Result of solving a `CoupledMPMProblem`. `component_matrices[name][step]` is the projection
matrix used for a component at a step; the vectors are empty if the problem was solved with
`save_matrices = false`.
"""
struct CoupledMPMSolution{T, U, M, S, O, R} <: AbstractProjectionSolution
    t::Vector{Int}
    u::U
    component_names::Vector{Symbol}
    component_matrices::M
    component_totals::Dict{Symbol, Vector{T}}
    substep_log::S
    observables::O
    event_log::Vector{Tuple{Int, Symbol}}
    rule_log::R
    retcode::Symbol
end

const StateDependentMPMSolution = CoupledMPMSolution
const HybridMPMSolution = CoupledMPMSolution

Base.getindex(sol::CoupledMPMSolution, name::Symbol) = sol.component_totals[name]

function Base.show(io::IO, sol::CoupledMPMSolution)
    print(io,
        "CoupledMPMSolution($(length(sol.component_names)) components, ",
        "$(length(sol.t)) time points, retcode=$(sol.retcode))")
end

# ---------------------------------------------------------------------------
# Solver helpers
# ---------------------------------------------------------------------------

function _materialize_component_matrix(comp::PopulationComponent, sys::PopulationSystem, day::Int, p)
    model = comp.model
    matrix = if model isa MatrixProjectionModel
        model.A
    elseif model isa AbstractMatrix
        model
    elseif applicable(model, sys, day, p)
        model(sys, day, p)
    else
        throw(ArgumentError(
            "component model $(typeof(model)) is neither a matrix nor callable as (system, day, p)"))
    end

    matrix isa MatrixProjectionModel && return _materialize_component_matrix(
        PopulationComponent(matrix, comp.population;
            stage_names = comp.stage_names,
            species = comp.species,
            type = comp.type,
            patch = comp.patch),
        sys,
        day,
        p)

    matrix isa AbstractMatrix || throw(ArgumentError(
        "component model must materialize to an AbstractMatrix or MatrixProjectionModel, got $(typeof(matrix))"))

    n = length(comp.population)
    size(matrix) == (n, n) || throw(DimensionMismatch(
        "component matrix has size $(size(matrix)); expected ($n, $n)"))
    return matrix
end

# One time step of a coupled matrix problem, on the shared engine.
_advance_one_day!(sys::PopulationSystem, day::Int, p, events, substeps, rules; kwargs...) =
    advance_coupled!(_materialize_component_matrix, sys, day, p, events, substeps, rules; kwargs...)

# ---------------------------------------------------------------------------
# Solve dispatch
# ---------------------------------------------------------------------------

function CommonSolve.solve(prob::CoupledMPMProblem, ::EigenAnalysis; kwargs...)
    throw(ArgumentError("EigenAnalysis is not defined for CoupledMPMProblem"))
end

"""
    solve(prob::CoupledMPMProblem, DirectIteration(); save_matrices = true)

Iterate a coupled problem one time step at a time. The [`CoupledMPMSolution`](@ref) holds the
state at every time point (`u`), the component totals, the substep and rule logs, the
observables and, for each component, the projection matrix used at each step
(`component_matrices`).

Keeping the matrices costs one dense matrix per component and time step. With
`save_matrices = false` they are not kept, and `component_matrices` holds an empty vector for
each component. Use it for long runs, or for components with many stages, when only the
trajectory is wanted.
"""
function CommonSolve.solve(prob::CoupledMPMProblem, ::DirectIteration = DirectIteration();
        save_matrices::Bool = true, kwargs...)
    r = solve_coupled(_materialize_component_matrix, prob.system, prob.tspan, prob.p,
        prob.substeps, prob.rules, prob.events, prob.observables, prob.normalize;
        save_matrices = save_matrices)
    return CoupledMPMSolution(r.t, r.u, r.component_names, r.component_matrices, r.component_totals,
        r.substep_log, r.observables, r.event_log, r.rule_log, r.retcode)
end
