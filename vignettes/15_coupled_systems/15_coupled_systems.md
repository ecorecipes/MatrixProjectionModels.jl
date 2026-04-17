# Coupled Population Systems
Simon Frost

## Overview

The coupled API lets multiple `PopulationComponent`s share state and
rules within a single `PopulationSystem`. This is the substrate for
hybrid matrix/PBDM workflows, multi-species interactions, and
patch-structured systems. This vignette covers component construction,
system introspection, state mutation, the priority-allocation pool,
observables, and the `CoupledMPMProblem` / `CoupledMPMSolution` types
together with their `Hybrid*` and `StateDependent*` aliases.

## Setup

``` julia
using MatrixProjectionModels
using LinearAlgebra
```

## Components

A `PopulationComponent` wraps a transition matrix (or an
`(sys, day, p) -> matrix` closure) plus a population vector and three
orthogonal labels: `species`, `type`, `patch`.

``` julia
A_wild = [0.6 2.0;
          0.4 0.5]
A_pred = [0.7 0.0;
          0.3 0.6]

wild = PopulationComponent(A_wild, [10.0, 5.0];
    stage_names = [:juv, :adult], species = :prey,    type = :wild,    patch = :north)
pred = PopulationComponent(A_pred, [2.0, 1.0];
    stage_names = [:juv, :adult], species = :predator, type = :wild,    patch = :north)

println(wild)
println(pred)
```

    PopulationComponent{Float64, Matrix{Float64}, Vector{Float64}}([0.6 2.0; 0.4 0.5], [10.0, 5.0], [:juv, :adult], :prey, :wild, :north)
    PopulationComponent{Float64, Matrix{Float64}, Vector{Float64}}([0.7 0.0; 0.3 0.6], [2.0, 1.0], [:juv, :adult], :predator, :wild, :north)

## Systems

`PopulationSystem` glues components together by name and carries shared
state.

``` julia
sys = PopulationSystem(:wild => wild, :pred => pred;
    state = (resource = 100.0, mortality_rate = 0.01))
println(sys)
```

    PopulationSystem(2 components: wild, pred, 2 state vars)

Total and per-component sizes:

``` julia
println("total_population   = ", total_population(sys))
println("component_total(wild) = ", component_total(sys, :wild))
println("component_totals     = ", component_totals(sys))
```

    total_population   = 18.0
    component_total(wild) = 15.0
    component_totals     = Dict(:pred => 3.0, :wild => 15.0)

State accessors:

``` julia
println("has_state(:resource)         = ", has_state(sys, :resource))
println("get_state(:resource)         = ", get_state(sys, :resource))
set_state!(sys, :resource, 95.0)
println("after set_state!, resource    = ", get_state(sys, :resource))
```

    has_state(:resource)         = true
    get_state(:resource)         = 100.0
    after set_state!, resource    = 95.0

Filters by species / type / patch return `(name, component)` pairs:

``` julia
println("by_species(:prey)     => ", [name for (name, _) in by_species(sys, :prey)])
println("by_type(:wild)        => ", [name for (name, _) in by_type(sys, :wild)])
println("by_patch(:north)      => ", [name for (name, _) in by_patch(sys, :north)])
```

    by_species(:prey)     => [:wild]
    by_type(:wild)        => [:wild, :pred]
    by_patch(:north)      => [:wild, :pred]

## Mutating populations

`inject!` adds individuals into a stage; `remove_fraction!` removes a
fraction of the population (whole component or a single stage).

``` julia
inject!(sys, :wild, 1, 3.0)
println("after inject!,    wild  = ", sys[:wild].population)
remove_fraction!(sys, :pred, 0.5)
println("after remove half, pred = ", sys[:pred].population)
```

    after inject!,    wild  = [13.0, 5.0]
    after remove half, pred = [1.0, 0.5]

## Priority allocation

`PriorityAllocationPool(supply, demands, labels)` allocates a finite
supply to ordered demands; later entries get whatever remains.

``` julia
pool = PriorityAllocationPool(10.0, [4.0, 5.0, 6.0], [:growth, :repro, :defence])
alloc = allocate(pool)
println("allocate            = ", alloc)
sdi   = supply_demand_index(pool)
println("supply_demand_index = ", round(sdi, digits=4))
allocs, stress = allocation_stress(pool)
println("allocations = ", allocs)
println("stress      = ", round.(stress, digits=4))
```

    allocate            = [4.0, 5.0, 1.0]
    supply_demand_index = 0.6667
    allocations = [4.0, 5.0, 1.0]
    stress      = [0.0, 0.0, 0.8333]

## Observables

An `Observable` is a `(sys, day, p) -> value` summary that is logged
each step.

``` julia
obs_total = Observable(:N_total, (sys, day, p) -> total_population(sys))
obs_pred  = Observable(:N_pred,  (sys, day, p) -> component_total(sys, :pred))
println("typeof(obs_total).name.name = ", typeof(obs_total).name.name)
```

    typeof(obs_total).name.name = Observable

## CoupledMPMProblem and aliases

``` julia
prob = CoupledMPMProblem(sys, (0, 5); observables = [obs_total, obs_pred])
println(prob)
```

    CoupledMPMProblem(2 components, substeps=0, rules=0, events=0, tspan=(0, 5))

`HybridMPMProblem` and `StateDependentMPMProblem` are exported aliases
for the same type — different names emphasising different intended
usage.

``` julia
println("HybridMPMProblem        === CoupledMPMProblem  = ",
        HybridMPMProblem        === CoupledMPMProblem)
println("StateDependentMPMProblem === CoupledMPMProblem = ",
        StateDependentMPMProblem === CoupledMPMProblem)
```

    HybridMPMProblem        === CoupledMPMProblem  = true
    StateDependentMPMProblem === CoupledMPMProblem = true

`remake` returns a new problem with overrides applied (for parameter
sweeps or re-running with different events):

``` julia
prob2 = remake(prob; tspan = (0, 3))
println("remade tspan = ", prob2.tspan)
```

    remade tspan = (0, 3)

## Solving and inspecting

``` julia
sol = solve(prob, DirectIteration())
println(sol)

println("retcode    = ", sol.retcode)
println("time steps = ", length(sol.t))
println("N_total[end] = ", round(sol.observables[:N_total][end], digits=4))
println("HybridMPMSolution        === CoupledMPMSolution = ",
        HybridMPMSolution        === CoupledMPMSolution)
println("StateDependentMPMSolution === CoupledMPMSolution = ",
        StateDependentMPMSolution === CoupledMPMSolution)
```

    CoupledMPMSolution(2 components, 6 time points, retcode=Success)
    retcode    = Success
    time steps = 6
    N_total[end] = 112.3489
    HybridMPMSolution        === CoupledMPMSolution = true
    StateDependentMPMSolution === CoupledMPMSolution = true

Per-component trajectory access:

``` julia
println("wild totals over time = ", round.(sol[:wild], digits=2))
```

    wild totals over time = [18.0, 25.5, 37.05, 53.5, 77.38, 111.87]

## Summary

- `PopulationComponent`, `PopulationSystem`, `Observable` are the core
  containers.
- `total_population`, `component_total`, `component_totals`,
  `by_species`, `by_type`, `by_patch` summarise.
- `inject!`, `remove_fraction!`, `get_state`, `set_state!`, `has_state`
  mutate / inspect.
- `PriorityAllocationPool` + `allocate` + `allocation_stress` +
  `supply_demand_index` provide resource-allocation building blocks.
- `CoupledMPMProblem` (alias `HybridMPMProblem`,
  `StateDependentMPMProblem`) and `CoupledMPMSolution` (alias
  `HybridMPMSolution`, `StateDependentMPMSolution`) drive the
  simulation; `remake` builds variants.
