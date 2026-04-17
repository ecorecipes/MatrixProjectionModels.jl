"""
Stateful discrete matrix population systems.

This layer complements `MPMProblem` for models that need multiple named
components, auxiliary mutable state, scheduled events, and post-step
transition rules without hand-written outer iteration loops.
"""

# ---------------------------------------------------------------------------
# Population components and systems
# ---------------------------------------------------------------------------

"""
    PopulationComponent(model, population; stage_names, species, type, patch)

A named stage-structured population component in a `PopulationSystem`.

- `model` may be a `MatrixProjectionModel`, an `AbstractMatrix`, or a callable
  `(system, day, p) -> AbstractMatrix`.
- `population` is the current stage vector.
- `species`, `type`, and `patch` carry lightweight grouping metadata.
"""
struct PopulationComponent{T<:Real, M, V<:AbstractVector{T}}
    model::M
    population::V
    stage_names::Vector{Symbol}
    species::Symbol
    type::Symbol
    patch::Symbol
end

function PopulationComponent(model, population::AbstractVector{<:Real};
        stage_names::AbstractVector{Symbol} = Symbol[],
        species::Symbol = :default,
        type::Symbol = :default,
        patch::Symbol = :default)
    pop = collect(float.(population))
    isempty(pop) && throw(ArgumentError("population must contain at least one stage"))
    inferred_stage_names = isempty(stage_names) ?
        _component_stage_names(model, length(pop)) :
        collect(stage_names)
    length(inferred_stage_names) == length(pop) || throw(DimensionMismatch(
        "stage_names length $(length(inferred_stage_names)) does not match population length $(length(pop))"))
    _validate_component_model(model, length(pop))
    return PopulationComponent{eltype(pop), typeof(model), typeof(pop)}(
        model,
        pop,
        inferred_stage_names,
        species,
        type,
        patch)
end

function _component_stage_names(model::MatrixProjectionModel, n::Int)
    size(model.A) == (n, n) || throw(DimensionMismatch(
        "model has size $(size(model.A)); expected ($n, $n)"))
    return copy(model.stage_names)
end

function _component_stage_names(model::AbstractMatrix, n::Int)
    size(model) == (n, n) || throw(DimensionMismatch(
        "model has size $(size(model)); expected ($n, $n)"))
    return _default_stage_names(n)
end

_component_stage_names(model, n::Int) = _default_stage_names(n)

function _validate_component_model(model::MatrixProjectionModel, n::Int)
    size(model.A) == (n, n) || throw(DimensionMismatch(
        "model has size $(size(model.A)); expected ($n, $n)"))
end

function _validate_component_model(model::AbstractMatrix, n::Int)
    size(model) == (n, n) || throw(DimensionMismatch(
        "model has size $(size(model)); expected ($n, $n)"))
end

_validate_component_model(model, n::Int) = nothing

function _convert_component(comp::PopulationComponent, ::Type{T}) where {T<:Real}
    PopulationComponent(comp.model, T.(comp.population);
        stage_names = comp.stage_names,
        species = comp.species,
        type = comp.type,
        patch = comp.patch)
end

function _copy_component(comp::PopulationComponent)
    PopulationComponent(comp.model, copy(comp.population);
        stage_names = comp.stage_names,
        species = comp.species,
        type = comp.type,
        patch = comp.patch)
end

"""
    PopulationSystem(components...; state = NamedTuple())

Ordered collection of named `PopulationComponent`s plus optional mutable
auxiliary state exposed to rules, events, and observables.
"""
struct PopulationSystem{T<:Real}
    order::Vector{Symbol}
    components::Dict{Symbol, PopulationComponent{T}}
    state::Dict{Symbol, Any}
end

function PopulationSystem(pairs::Pair{Symbol}...; state = NamedTuple())
    isempty(pairs) && throw(ArgumentError("PopulationSystem must have at least one component"))

    order = Symbol[]
    raw_components = Pair{Symbol, Any}[]
    for (name, value) in pairs
        name in order && throw(ArgumentError("duplicate component name: $name"))
        value isa PopulationComponent || throw(ArgumentError(
            "expected PopulationComponent for :$name, got $(typeof(value))"))
        push!(order, name)
        push!(raw_components, name => value)
    end

    T = promote_type(map(p -> eltype(p.second.population), raw_components)...)
    components = Dict{Symbol, PopulationComponent{T}}()
    for (name, comp) in raw_components
        components[name] = _convert_component(comp, T)
    end

    return PopulationSystem{T}(copy(order), components, _coerce_state(state))
end

function _coerce_state(state)
    if state isa NamedTuple
        return Dict{Symbol, Any}(k => deepcopy(v) for (k, v) in pairs(state))
    elseif state isa AbstractDict
        return Dict{Symbol, Any}(Symbol(k) => deepcopy(v) for (k, v) in pairs(state))
    elseif state isa AbstractVector{<:Pair}
        return Dict{Symbol, Any}(Symbol(k) => deepcopy(v) for (k, v) in state)
    else
        throw(ArgumentError("state must be a NamedTuple, AbstractDict, or vector of pairs"))
    end
end

Base.getindex(sys::PopulationSystem, name::Symbol) = sys.components[name]
Base.haskey(sys::PopulationSystem, name::Symbol) = haskey(sys.components, name)
Base.keys(sys::PopulationSystem) = sys.order
Base.length(sys::PopulationSystem) = length(sys.order)
Base.pairs(sys::PopulationSystem) = ((name, sys.components[name]) for name in sys.order)

function Base.show(io::IO, sys::PopulationSystem)
    state_suffix = isempty(sys.state) ? "" : ", $(length(sys.state)) state vars"
    print(io, "PopulationSystem($(length(sys)) components: $(join(sys.order, ", "))$state_suffix)")
end

get_state(sys::PopulationSystem, name::Symbol) = sys.state[name]
set_state!(sys::PopulationSystem, name::Symbol, value) = (sys.state[name] = value)
has_state(sys::PopulationSystem, name::Symbol) = haskey(sys.state, name)

function _copy_system(sys::PopulationSystem)
    return PopulationSystem(
        (name => _copy_component(sys[name]) for name in sys.order)...;
        state = deepcopy(sys.state))
end

total_population(sys::PopulationSystem) = sum(sum(comp.population) for (_, comp) in pairs(sys))
component_total(sys::PopulationSystem, name::Symbol) = sum(sys[name].population)
component_totals(sys::PopulationSystem{T}) where {T} =
    Dict{Symbol, T}(name => sum(comp.population) for (name, comp) in pairs(sys))

by_species(sys::PopulationSystem, species::Symbol) =
    [(name, comp) for (name, comp) in pairs(sys) if comp.species == species]

by_type(sys::PopulationSystem, type::Symbol) =
    [(name, comp) for (name, comp) in pairs(sys) if comp.type == type]

by_patch(sys::PopulationSystem, patch::Symbol) =
    [(name, comp) for (name, comp) in pairs(sys) if comp.patch == patch]

function _validate_stage(comp::PopulationComponent, stage::Int)
    1 <= stage <= length(comp.population) || throw(BoundsError(
        "stage $stage out of range for component with $(length(comp.population)) stages"))
end

function inject!(sys::PopulationSystem, name::Symbol, stage::Int, amount::Real)
    amount >= 0 || throw(ArgumentError("injected amount must be nonnegative"))
    comp = sys[name]
    _validate_stage(comp, stage)
    comp.population[stage] += amount
    return nothing
end

inject!(sys::PopulationSystem, name::Symbol, amount::Real) = inject!(sys, name, 1, amount)

function remove_fraction!(sys::PopulationSystem, name::Symbol, fraction::Real)
    _validate_fraction("fraction", fraction)
    comp = sys[name]
    comp.population .*= (one(eltype(comp.population)) - eltype(comp.population)(fraction))
    return nothing
end

function remove_fraction!(sys::PopulationSystem, name::Symbol, stage::Int, fraction::Real)
    _validate_fraction("fraction", fraction)
    comp = sys[name]
    _validate_stage(comp, stage)
    comp.population[stage] *= (one(eltype(comp.population)) - eltype(comp.population)(fraction))
    return nothing
end

# ---------------------------------------------------------------------------
# Hybrid allocation helpers and ordered daily substeps
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

"""
    AbstractDailySubstep

Ordered callbacks executed after scheduled events and before matrix
materialization. Substeps mutate auxiliary state and/or component populations to
support hybrid daily workflows.
"""
abstract type AbstractDailySubstep end

function apply_substep! end

"""
    StateUpdateSubstep(name, state_key, update_fn)

Set `system[state_key] = update_fn(system, day, p)` and record the new value.
"""
struct StateUpdateSubstep{F} <: AbstractDailySubstep
    name::Symbol
    state_key::Symbol
    update_fn::F
end

function apply_substep!(step::StateUpdateSubstep, sys::PopulationSystem, day::Int, p)
    value = deepcopy(step.update_fn(sys, day, p))
    set_state!(sys, step.state_key, value)
    return NamedTuple{(step.state_key,)}((deepcopy(value),))
end

"""
    CustomSubstep(name, apply_fn)

Arbitrary ordered pre-step callback. `apply_fn(system, day, p)` may mutate the
system and should return either a `NamedTuple` of metrics or `nothing`.
"""
struct CustomSubstep{F} <: AbstractDailySubstep
    name::Symbol
    apply_fn::F
end

function apply_substep!(step::CustomSubstep, sys::PopulationSystem, day::Int, p)
    result = step.apply_fn(sys, day, p)
    result === nothing && return NamedTuple()
    result isa NamedTuple || throw(ArgumentError(
        "custom substeps must return a NamedTuple or nothing, got $(typeof(result))"))
    return deepcopy(result)
end

_substep_name(step::AbstractDailySubstep) = getfield(step, :name)

# ---------------------------------------------------------------------------
# Fluxes, rules, events, and observables
# ---------------------------------------------------------------------------

"""
    StageFlux(component, stage, amount)

Signed stage-targeted population change applied after rule evaluation.
Positive values add individuals; negative values remove them.
"""
struct StageFlux{T<:Real}
    component::Symbol
    stage::Int
    amount::T

    function StageFlux{T}(component::Symbol, stage::Int, amount::T) where {T<:Real}
        stage >= 1 || throw(ArgumentError("stage must be >= 1"))
        isfinite(amount) || throw(ArgumentError("amount must be finite"))
        new{T}(component, stage, amount)
    end
end

StageFlux(component::Symbol, stage::Int, amount::Real) =
    StageFlux{typeof(float(amount))}(component, stage, float(amount))

"""
    StageTransfer(source, source_stage, target, target_stage, amount)

Stage-targeted transfer applied after rule evaluation with mass-balance checks.
"""
struct StageTransfer{T<:Real}
    source::Symbol
    source_stage::Int
    target::Symbol
    target_stage::Int
    amount::T

    function StageTransfer{T}(source::Symbol, source_stage::Int,
            target::Symbol, target_stage::Int, amount::T) where {T<:Real}
        source_stage >= 1 || throw(ArgumentError("source_stage must be >= 1"))
        target_stage >= 1 || throw(ArgumentError("target_stage must be >= 1"))
        isfinite(amount) || throw(ArgumentError("amount must be finite"))
        amount >= zero(T) || throw(ArgumentError("transfer amount must be nonnegative"))
        new{T}(source, source_stage, target, target_stage, amount)
    end
end

StageTransfer(source::Symbol, source_stage::Int, target::Symbol, target_stage::Int, amount::Real) =
    StageTransfer{typeof(float(amount))}(source, source_stage, target, target_stage, float(amount))

"""
    RuleEffect(; fluxes = StageFlux[], transfers = StageTransfer[], metrics = NamedTuple())

Normalized output of a transition rule.
"""
struct RuleEffect{T<:Real, M<:NamedTuple}
    fluxes::Vector{StageFlux{T}}
    transfers::Vector{StageTransfer{T}}
    metrics::M
end

function RuleEffect(; fluxes::AbstractVector = StageFlux[],
        transfers::AbstractVector = StageTransfer[],
        metrics = NamedTuple())
    metrics isa NamedTuple || throw(ArgumentError("metrics must be a NamedTuple"))

    amount_types = Type[]
    for flux in fluxes
        flux isa StageFlux || throw(ArgumentError("fluxes must contain StageFlux values"))
        push!(amount_types, typeof(flux.amount))
    end
    for transfer in transfers
        transfer isa StageTransfer || throw(ArgumentError("transfers must contain StageTransfer values"))
        push!(amount_types, typeof(transfer.amount))
    end

    T = isempty(amount_types) ? Float64 : promote_type(amount_types...)
    flux_vec = StageFlux{T}[StageFlux(flux.component, flux.stage, T(flux.amount)) for flux in fluxes]
    transfer_vec = StageTransfer{T}[StageTransfer(transfer.source, transfer.source_stage,
        transfer.target, transfer.target_stage, T(transfer.amount)) for transfer in transfers]
    return RuleEffect{T, typeof(metrics)}(flux_vec, transfer_vec, metrics)
end

function _normalize_rule_effect(result)
    result === nothing && return RuleEffect()
    result isa RuleEffect && return result
    result isa StageFlux && return RuleEffect(fluxes = [result])
    result isa StageTransfer && return RuleEffect(transfers = [result])

    if result isa AbstractVector
        all(x -> x isa StageFlux, result) && return RuleEffect(fluxes = result)
        all(x -> x isa StageTransfer, result) && return RuleEffect(transfers = result)
        throw(ArgumentError("vectors returned by rules must contain only StageFlux or only StageTransfer values"))
    end

    if result isa NamedTuple
        names = propertynames(result)
        if any(name -> name in (:fluxes, :transfers, :metrics), names)
            fluxes = hasproperty(result, :fluxes) ? getproperty(result, :fluxes) : StageFlux[]
            transfers = hasproperty(result, :transfers) ? getproperty(result, :transfers) : StageTransfer[]
            metrics = hasproperty(result, :metrics) ? getproperty(result, :metrics) : NamedTuple()
            return RuleEffect(; fluxes = fluxes, transfers = transfers, metrics = metrics)
        end
        return RuleEffect(metrics = result)
    end

    throw(ArgumentError(
        "rule outputs must be RuleEffect, StageFlux, StageTransfer, vectors of them, or NamedTuples"))
end

abstract type AbstractTransitionRule end
abstract type AbstractScheduledEvent end

function apply_rule end
function apply_event! end

_validate_fraction(label::AbstractString, value::Real) =
    (isfinite(value) && 0 <= value <= 1) ||
    throw(ArgumentError("$label must be finite and lie in [0, 1], got $value"))

function _validate_nonnegative(label::AbstractString, value::Real)
    isfinite(value) || throw(ArgumentError("$label must be finite, got $value"))
    value >= 0 || throw(ArgumentError("$label must be nonnegative, got $value"))
    return value
end

"""
    TransferRule(source, target, fraction_fn; source_stage, target_stage, name)

Transfer a fraction of one source stage into a target stage after stepping.
"""
struct TransferRule{F} <: AbstractTransitionRule
    name::Symbol
    source::Symbol
    source_stage::Int
    target::Symbol
    target_stage::Int
    fraction_fn::F
end

function TransferRule(source::Symbol, target::Symbol, fraction_fn;
        source_stage::Int = 1,
        target_stage::Int = 1,
        name::Symbol = Symbol(source, "_to_", target))
    return TransferRule(name, source, source_stage, target, target_stage, fraction_fn)
end

function apply_rule(rule::TransferRule, sys::PopulationSystem, day::Int, p)
    comp = sys[rule.source]
    _validate_stage(comp, rule.source_stage)
    _validate_stage(sys[rule.target], rule.target_stage)
    frac = rule.fraction_fn(sys, day, p)
    _validate_fraction("transfer fraction", frac)
    amount = comp.population[rule.source_stage] * frac
    return RuleEffect(
        transfers = [StageTransfer(rule.source, rule.source_stage, rule.target, rule.target_stage, amount)],
        metrics = (transferred = amount,))
end

"""
    ReproductionRule(target, reproduction_fn; stage, name)

Inject offspring into a target stage after stepping.
"""
struct ReproductionRule{F} <: AbstractTransitionRule
    name::Symbol
    target::Symbol
    stage::Int
    reproduction_fn::F
end

function ReproductionRule(target::Symbol, reproduction_fn;
        stage::Int = 1,
        name::Symbol = Symbol("reproduce_", target))
    stage >= 1 || throw(ArgumentError("stage must be >= 1"))
    return ReproductionRule(name, target, stage, reproduction_fn)
end

function apply_rule(rule::ReproductionRule, sys::PopulationSystem, day::Int, p)
    _validate_stage(sys[rule.target], rule.stage)
    offspring = rule.reproduction_fn(sys, day, p)
    _validate_nonnegative("offspring", offspring)
    return RuleEffect(
        fluxes = [StageFlux(rule.target, rule.stage, offspring)],
        metrics = (offspring = offspring,))
end

"""
    MortalityRule(target, mortality_fn; stage, name)

Apply additional post-step mortality to one stage or all stages (`stage = 0`).
"""
struct MortalityRule{F} <: AbstractTransitionRule
    name::Symbol
    target::Symbol
    stage::Int
    mortality_fn::F
end

function MortalityRule(target::Symbol, mortality_fn;
        stage::Int = 0,
        name::Symbol = Symbol("mortality_", target))
    stage >= 0 || throw(ArgumentError("stage must be >= 0"))
    return MortalityRule(name, target, stage, mortality_fn)
end

function apply_rule(rule::MortalityRule, sys::PopulationSystem, day::Int, p)
    comp = sys[rule.target]
    frac = rule.mortality_fn(sys, day, p)
    _validate_fraction("mortality fraction", frac)

    fluxes = StageFlux[]
    if rule.stage == 0
        for stage in eachindex(comp.population)
            push!(fluxes, StageFlux(rule.target, stage, -comp.population[stage] * frac))
        end
    else
        _validate_stage(comp, rule.stage)
        push!(fluxes, StageFlux(rule.target, rule.stage, -comp.population[rule.stage] * frac))
    end

    removed = -sum(flux.amount for flux in fluxes)
    return RuleEffect(fluxes = fluxes, metrics = (removed = removed, mortality = frac))
end

"""
    CustomRule(name, apply_fn)

Escape hatch for user-defined rule logic. `apply_fn(system, day, p)` must return
one of the accepted `RuleEffect`-compatible outputs.
"""
struct CustomRule{F} <: AbstractTransitionRule
    name::Symbol
    apply_fn::F
end

apply_rule(rule::CustomRule, sys::PopulationSystem, day::Int, p) =
    _normalize_rule_effect(rule.apply_fn(sys, day, p))

"""
    PulseRelease(target, amount, interval; stage_idx, start_day, end_day, name)

Inject a fixed amount on a regular day interval.
"""
struct PulseRelease{T<:Real} <: AbstractScheduledEvent
    name::Symbol
    target::Symbol
    stage_idx::Int
    amount::T
    interval::Int
    start_day::Int
    end_day::Int
end

function PulseRelease(target::Symbol, amount::Real, interval::Int;
        stage_idx::Int = 1,
        start_day::Int = 0,
        end_day::Int = typemax(Int),
        name::Symbol = Symbol("pulse_", target))
    interval > 0 || throw(ArgumentError("interval must be positive"))
    stage_idx >= 1 || throw(ArgumentError("stage_idx must be >= 1"))
    _validate_nonnegative("release amount", amount)
    return PulseRelease(name, target, stage_idx, float(amount), interval, start_day, end_day)
end

function apply_event!(ev::PulseRelease, sys::PopulationSystem, day::Int, p)
    day < ev.start_day && return false
    day > ev.end_day && return false
    (day - ev.start_day) % ev.interval == 0 || return false
    inject!(sys, ev.target, ev.stage_idx, ev.amount)
    return true
end

"""
    SingleDayRelease(target, amount, day; stage_idx, name)

Inject a fixed amount on exactly one day.
"""
struct SingleDayRelease{T<:Real} <: AbstractScheduledEvent
    name::Symbol
    target::Symbol
    stage_idx::Int
    amount::T
    day::Int
end

function SingleDayRelease(target::Symbol, amount::Real, day::Int;
        stage_idx::Int = 1,
        name::Symbol = Symbol("release_", target))
    stage_idx >= 1 || throw(ArgumentError("stage_idx must be >= 1"))
    _validate_nonnegative("release amount", amount)
    return SingleDayRelease(name, target, stage_idx, float(amount), day)
end

function apply_event!(ev::SingleDayRelease, sys::PopulationSystem, day::Int, p)
    day == ev.day || return false
    inject!(sys, ev.target, ev.stage_idx, ev.amount)
    return true
end

"""
    SprayEvent(targets, kill_fractions, days; name)

Apply multiplicative mortality to named components on specified days.
"""
struct SprayEvent <: AbstractScheduledEvent
    name::Symbol
    targets::Vector{Symbol}
    kill_fractions::Vector{Float64}
    days::Vector{Int}
end

function SprayEvent(targets::AbstractVector{Symbol},
        kill_fractions::AbstractVector{<:Real},
        days::AbstractVector{Int};
        name::Symbol = :spray)
    length(targets) == length(kill_fractions) || throw(DimensionMismatch(
        "targets length $(length(targets)) does not match kill_fractions length $(length(kill_fractions))"))
    foreach(frac -> _validate_fraction("kill fraction", frac), kill_fractions)
    return SprayEvent(name, collect(targets), Float64.(kill_fractions), collect(days))
end

function apply_event!(ev::SprayEvent, sys::PopulationSystem, day::Int, p)
    day in ev.days || return false
    for (target, frac) in zip(ev.targets, ev.kill_fractions)
        haskey(sys, target) || throw(ArgumentError("unknown spray target :$target"))
        remove_fraction!(sys, target, frac)
    end
    return true
end

"""
    CustomEvent(name, apply_fn)

Custom scheduled event. `apply_fn(system, day, p)` must mutate the system if it
fires and return `true`; otherwise return `false`.
"""
struct CustomEvent{F} <: AbstractScheduledEvent
    name::Symbol
    apply_fn::F
end

function apply_event!(ev::CustomEvent, sys::PopulationSystem, day::Int, p)
    result = ev.apply_fn(sys, day, p)
    result isa Bool || throw(ArgumentError("custom events must return Bool, got $(typeof(result))"))
    return result
end

"""
    Observable(name, fn)

Record a scalar observable at each stored time point.
"""
struct Observable{F}
    name::Symbol
    fn::F
end

_rule_name(rule::AbstractTransitionRule) = getfield(rule, :name)
_event_name(event::AbstractScheduledEvent) = getfield(event, :name)

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

function _ensure_unique_names(items, label::AbstractString)
    names = Symbol[]
    for item in items
        hasfield(typeof(item), :name) || continue
        name = getfield(item, :name)
        name in names && throw(ArgumentError("duplicate $label name: $name"))
        push!(names, name)
    end
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

Result of solving a `CoupledMPMProblem`.
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

function _snapshot(sys::PopulationSystem{T}) where {T}
    return Dict{Symbol, Vector{T}}(name => copy(sys[name].population) for name in sys.order)
end

function _apply_rule_effects!(sys::PopulationSystem{T}, effects::AbstractVector{<:RuleEffect}) where {T}
    deltas = Dict{Symbol, Vector{T}}(name => zeros(T, length(sys[name].population)) for name in sys.order)

    for effect in effects
        for transfer in effect.transfers
            haskey(sys, transfer.source) || throw(ArgumentError("unknown transfer source :$(transfer.source)"))
            haskey(sys, transfer.target) || throw(ArgumentError("unknown transfer target :$(transfer.target)"))
            _validate_stage(sys[transfer.source], transfer.source_stage)
            _validate_stage(sys[transfer.target], transfer.target_stage)
            amount = T(transfer.amount)
            deltas[transfer.source][transfer.source_stage] -= amount
            deltas[transfer.target][transfer.target_stage] += amount
        end

        for flux in effect.fluxes
            haskey(sys, flux.component) || throw(ArgumentError("unknown flux component :$(flux.component)"))
            _validate_stage(sys[flux.component], flux.stage)
            deltas[flux.component][flux.stage] += T(flux.amount)
        end
    end

    for name in sys.order
        new_population = sys[name].population .+ deltas[name]
        any(value -> value < zero(T), new_population) && throw(ArgumentError(
            "rule effects would make component :$name negative"))
    end

    for name in sys.order
        sys[name].population .+= deltas[name]
    end

    return nothing
end

function _normalize_components!(sys::PopulationSystem)
    for name in sys.order
        pop = sys[name].population
        total = sum(pop)
        total > zero(total) || continue
        pop ./= total
    end
    return nothing
end

function _advance_one_day!(sys::PopulationSystem, day::Int, p, events, substeps, rules;
        normalize::Bool = false,
        event_log = nothing)
    for event in events
        fired = apply_event!(event, sys, day, p)
        if fired && event_log !== nothing
            push!(event_log, (day, _event_name(event)))
        end
    end

    substep_metrics = NamedTuple[]
    for substep in substeps
        metrics = apply_substep!(substep, sys, day, p)
        metrics isa NamedTuple || throw(ArgumentError(
            "substeps must return NamedTuple metrics, got $(typeof(metrics))"))
        push!(substep_metrics, deepcopy(metrics))
    end

    current = _snapshot(sys)
    matrices = Dict{Symbol, Any}()
    stepped = Dict{Symbol, Any}()

    for name in sys.order
        comp = sys[name]
        A = _materialize_component_matrix(comp, sys, day, p)
        matrices[name] = Matrix(A)
        stepped[name] = Vector(A * current[name])
    end

    for name in sys.order
        sys[name].population .= stepped[name]
    end

    effects = RuleEffect[apply_rule(rule, sys, day, p) for rule in rules]
    _apply_rule_effects!(sys, effects)
    normalize && _normalize_components!(sys)

    return matrices, substep_metrics, effects
end

# ---------------------------------------------------------------------------
# Solve dispatch
# ---------------------------------------------------------------------------

function CommonSolve.solve(prob::CoupledMPMProblem, ::EigenAnalysis; kwargs...)
    throw(ArgumentError("EigenAnalysis is not defined for CoupledMPMProblem"))
end

function CommonSolve.solve(prob::CoupledMPMProblem, ::DirectIteration = DirectIteration(); kwargs...)
    sys = _copy_system(prob.system)
    t0, tf = prob.tspan
    n_steps = tf - t0
    ts = collect(t0:tf)

    snapshots = Vector{Dict{Symbol, Vector{eltype(first(values(sys.components)).population)}}}(undef, n_steps + 1)
    snapshots[1] = _snapshot(sys)

    component_totals_store = Dict{Symbol, Vector{eltype(first(values(sys.components)).population)}}(
        name => Vector{eltype(first(values(sys.components)).population)}(undef, n_steps + 1)
        for name in sys.order)
    for name in sys.order
        component_totals_store[name][1] = sum(sys[name].population)
    end

    observable_store = Dict{Symbol, Vector{Any}}(
        obs.name => Vector{Any}(undef, n_steps + 1) for obs in prob.observables)
    for obs in prob.observables
        observable_store[obs.name][1] = obs.fn(sys, t0, prob.p)
    end

    component_matrices = Dict{Symbol, Vector{Any}}(
        name => Vector{Any}(undef, n_steps) for name in sys.order)
    substep_log = Dict{Symbol, Vector{Any}}(
        _substep_name(substep) => Vector{Any}(undef, n_steps) for substep in prob.substeps)
    rule_log = Dict{Symbol, Vector{Any}}(
        _rule_name(rule) => Vector{Any}(undef, n_steps) for rule in prob.rules)
    event_log = Tuple{Int, Symbol}[]

    for step in 1:n_steps
        day = t0 + step - 1
        matrices, substep_metrics, effects = _advance_one_day!(
            sys,
            day,
            prob.p,
            prob.events,
            prob.substeps,
            prob.rules;
            normalize = prob.normalize,
            event_log = event_log)

        for name in sys.order
            component_matrices[name][step] = matrices[name]
            component_totals_store[name][step + 1] = sum(sys[name].population)
        end

        for (substep, metrics) in zip(prob.substeps, substep_metrics)
            substep_log[_substep_name(substep)][step] = metrics
        end

        for (rule, effect) in zip(prob.rules, effects)
            rule_log[_rule_name(rule)][step] = effect.metrics
        end

        snapshots[step + 1] = _snapshot(sys)
        for obs in prob.observables
            observable_store[obs.name][step + 1] = obs.fn(sys, ts[step + 1], prob.p)
        end
    end

    return CoupledMPMSolution(
        ts,
        snapshots,
        copy(sys.order),
        component_matrices,
        component_totals_store,
        substep_log,
        observable_store,
        event_log,
        rule_log,
        :Success)
end
