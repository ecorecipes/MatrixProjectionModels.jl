# Coupled Population Systems

## Overview

The coupled API lets multiple `PopulationComponent`s share state and rules
within a single `PopulationSystem`. This is the substrate for hybrid
matrix/PBDM workflows, multi-species interactions, and patch-structured
systems. This vignette covers component construction, system
introspection, state mutation, the priority-allocation pool, observables,
and the `CoupledMPMProblem` / `CoupledMPMSolution` types together with
their `Hybrid*` and `StateDependent*` aliases.

## Setup

```@example mpm
using MatrixProjectionModels
using LinearAlgebra
```

## Components

A `PopulationComponent` wraps a transition matrix (or an `(sys, day, p) ->
matrix` closure) plus a population vector and three orthogonal labels:
`species`, `type`, `patch`.

```@example mpm
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

## Systems

`PopulationSystem` glues components together by name and carries shared
state.

```@example mpm
sys = PopulationSystem(:wild => wild, :pred => pred;
    state = (resource = 100.0, mortality_rate = 0.01))
println(sys)
```

Total and per-component sizes:

```@example mpm
println("total_population   = ", total_population(sys))
println("component_total(wild) = ", component_total(sys, :wild))
println("component_totals     = ", component_totals(sys))
```

State accessors:

```@example mpm
println("has_state(:resource)         = ", has_state(sys, :resource))
println("get_state(:resource)         = ", get_state(sys, :resource))
set_state!(sys, :resource, 95.0)
println("after set_state!, resource    = ", get_state(sys, :resource))
```

Filters by species / type / patch return `(name, component)` pairs:

```@example mpm
println("by_species(:prey)     => ", [name for (name, _) in by_species(sys, :prey)])
println("by_type(:wild)        => ", [name for (name, _) in by_type(sys, :wild)])
println("by_patch(:north)      => ", [name for (name, _) in by_patch(sys, :north)])
```

## Mutating populations

`inject!` adds individuals into a stage; `remove_fraction!` removes a
fraction of the population (whole component or a single stage).

```@example mpm
inject!(sys, :wild, 1, 3.0)
println("after inject!,    wild  = ", sys[:wild].population)
remove_fraction!(sys, :pred, 0.5)
println("after remove half, pred = ", sys[:pred].population)
```

## Priority allocation

`PriorityAllocationPool(supply, demands, labels)` allocates a finite supply
to ordered demands; later entries get whatever remains.

```@example mpm
pool = PriorityAllocationPool(10.0, [4.0, 5.0, 6.0], [:growth, :repro, :defence])
alloc = allocate(pool)
println("allocate            = ", alloc)
sdi   = supply_demand_index(pool)
println("supply_demand_index = ", round(sdi, digits=4))
allocs, stress = allocation_stress(pool)
println("allocations = ", allocs)
println("stress      = ", round.(stress, digits=4))
```

## Observables

An `Observable` is a `(sys, day, p) -> value` summary that is logged each
step.

```@example mpm
obs_total = Observable(:N_total, (sys, day, p) -> total_population(sys))
obs_pred  = Observable(:N_pred,  (sys, day, p) -> component_total(sys, :pred))
println("typeof(obs_total).name.name = ", typeof(obs_total).name.name)
```

## CoupledMPMProblem and aliases

```@example mpm
prob = CoupledMPMProblem(sys, (0, 5); observables = [obs_total, obs_pred])
println(prob)
```

`HybridMPMProblem` and `StateDependentMPMProblem` are exported aliases for
the same type — different names emphasising different intended usage.

```@example mpm
println("HybridMPMProblem        === CoupledMPMProblem  = ",
        HybridMPMProblem        === CoupledMPMProblem)
println("StateDependentMPMProblem === CoupledMPMProblem = ",
        StateDependentMPMProblem === CoupledMPMProblem)
```

`remake` returns a new problem with overrides applied (for parameter sweeps
or re-running with different events):

```@example mpm
prob2 = remake(prob; tspan = (0, 3))
println("remade tspan = ", prob2.tspan)
```

## Solving and inspecting

```@example mpm
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

Per-component trajectory access:

```@example mpm
println("wild totals over time = ", round.(sol[:wild], digits=2))
```

## Summary

- `PopulationComponent`, `PopulationSystem`, `Observable` are the core
  containers.
- `total_population`, `component_total`, `component_totals`, `by_species`,
  `by_type`, `by_patch` summarise.
- `inject!`, `remove_fraction!`, `get_state`, `set_state!`, `has_state`
  mutate / inspect.
- `PriorityAllocationPool` + `allocate` + `allocation_stress` +
  `supply_demand_index` provide resource-allocation building blocks.
- `CoupledMPMProblem` (alias `HybridMPMProblem`, `StateDependentMPMProblem`)
  and `CoupledMPMSolution` (alias `HybridMPMSolution`,
  `StateDependentMPMSolution`) drive the simulation; `remake` builds
  variants.
