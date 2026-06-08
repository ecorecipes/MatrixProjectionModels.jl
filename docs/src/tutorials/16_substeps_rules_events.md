# Substeps, Rules, and Scheduled Events

## Overview

A `CoupledMPMProblem` advances each day in three phases:

1. **Substeps** (`AbstractDailySubstep`) update shared `state` *before*
   matrix multiplication.
2. **Matrix step** projects each component using its (possibly
   state-dependent) transition matrix.
3. **Rules** (`AbstractTransitionRule`) emit signed `StageFlux` /
   `StageTransfer` effects; **events** (`AbstractScheduledEvent`) mutate
   the system on specific days.

This vignette walks through every concrete subtype.

## Setup

```@example mpm
using MatrixProjectionModels
using LinearAlgebra
```

```@example mpm
A = [0.6 0.0;
     0.4 0.5]
comp = PopulationComponent(A, [10.0, 5.0]; stage_names = [:juv, :adult])
sys  = PopulationSystem(:pop => comp; state = (resource = 100.0, log = 0.0))
```

## Substeps

### `StateUpdateSubstep`

Updates a single named state variable.

```@example mpm
sub_resource = StateUpdateSubstep(:resource_decay, :resource,
    (sys, day, p) -> get_state(sys, :resource) - 1.0)
println("typeof(sub_resource).name.name = ", typeof(sub_resource).name.name)
println(sub_resource isa AbstractDailySubstep)
```

### `CustomSubstep`

Arbitrary mutation; the return value can be anything.

```@example mpm
sub_log = CustomSubstep(:log_step, (sys, day, p) -> begin
    set_state!(sys, :log, get_state(sys, :log) + 1.0)
    return (logged = day,)
end)
println("typeof(sub_log).name.name = ", typeof(sub_log).name.name)
```

`apply_substep!` is the dispatch hook; the solver calls it once per
substep per day.

```@example mpm
println("apply_substep! method count = ", length(methods(apply_substep!)))
```

## Rule effect primitives

`StageFlux` (signed scalar) and `StageTransfer` (mass-balance preserving
flow) are the building blocks; `RuleEffect` collects both plus a metrics
`NamedTuple`.

```@example mpm
flx  = StageFlux(:pop, 2, 1.5)
tfr  = StageTransfer(:pop, 1, :pop, 2, 0.4)
eff  = RuleEffect(fluxes = [flx], transfers = [tfr], metrics = (note = "demo",))
println("typeof(flx).name.name = ", typeof(flx).name.name)
println("typeof(tfr).name.name = ", typeof(tfr).name.name)
println("typeof(eff).name.name = ", typeof(eff).name.name)
```

## Transition rules

All four concrete rule types are `<: AbstractTransitionRule` and have the
shared `apply_rule(rule, sys, day, p)` interface.

```@example mpm
println("apply_rule method count = ", length(methods(apply_rule)))
```

### `TransferRule`

```@example mpm
tr = TransferRule(:pop, :pop, (sys, day, p) -> 0.1;
    source_stage = 1, target_stage = 2, name = :graduate)
println("typeof(tr).name.name = ", typeof(tr).name.name)
println(tr isa AbstractTransitionRule)
```

### `ReproductionRule`

```@example mpm
rr = ReproductionRule(:pop, (sys, day, p) -> 2.0 * sys[:pop].population[2];
    stage = 1, name = :spawn)
println("typeof(rr).name.name = ", typeof(rr).name.name)
```

### `MortalityRule`

```@example mpm
mr = MortalityRule(:pop, (sys, day, p) -> 0.05; stage = 0, name = :background)
println("typeof(mr).name.name = ", typeof(mr).name.name)
```

### `CustomRule`

```@example mpm
cr = CustomRule(:noop, (sys, day, p) -> RuleEffect())
println("typeof(cr).name.name = ", typeof(cr).name.name)
```

## Scheduled events

All four concrete event types are `<: AbstractScheduledEvent` and share
`apply_event!(event, sys, day, p)`.

```@example mpm
println("apply_event! method count = ", length(methods(apply_event!)))
```

### `PulseRelease`

Fixed amount on a regular interval.

```@example mpm
pulse = PulseRelease(:pop, 5.0, 2; stage_idx = 1, start_day = 0, end_day = 10,
    name = :weekly_pulse)
println("typeof(pulse).name.name = ", typeof(pulse).name.name)
println(pulse isa AbstractScheduledEvent)
```

### `SingleDayRelease`

```@example mpm
sdr = SingleDayRelease(:pop, 50.0, 3; stage_idx = 1, name = :one_off)
println("typeof(sdr).name.name = ", typeof(sdr).name.name)
```

### `SprayEvent`

Multiplicative mortality on listed days.

```@example mpm
spray = SprayEvent([:pop], [0.3], [4]; name = :spray)
println("typeof(spray).name.name = ", typeof(spray).name.name)
```

### `CustomEvent`

```@example mpm
ce = CustomEvent(:custom_pulse, (sys, day, p) -> begin
    if day == 1
        inject!(sys, :pop, 1, 2.0)
        return true
    else
        return false
    end
end)
println("typeof(ce).name.name = ", typeof(ce).name.name)
```

## End-to-end: assembling and solving

```@example mpm
prob = CoupledMPMProblem(sys, (0, 5);
    substeps    = [sub_resource, sub_log],
    rules       = [tr, rr, mr, cr],
    events      = [pulse, sdr, spray, ce],
    observables = [Observable(:N, (s, d, p) -> total_population(s)),
                   Observable(:R, (s, d, p) -> get_state(s, :resource))])

sol = solve(prob, DirectIteration())
println("retcode = ", sol.retcode)
println("days    = ", sol.t)
println("N       = ", round.(sol.observables[:N], digits=3))
println("R       = ", round.(sol.observables[:R], digits=3))
println("triggered events    = ", length(sol.event_log))
println("rule-log keys       = ", sort(collect(keys(sol.rule_log))))
```

## Summary

- Substeps: `AbstractDailySubstep`, `apply_substep!`, `StateUpdateSubstep`,
  `CustomSubstep`.
- Rule plumbing: `StageFlux`, `StageTransfer`, `RuleEffect`.
- Rules: `AbstractTransitionRule`, `apply_rule`, `TransferRule`,
  `ReproductionRule`, `MortalityRule`, `CustomRule`.
- Events: `AbstractScheduledEvent`, `apply_event!`, `PulseRelease`,
  `SingleDayRelease`, `SprayEvent`, `CustomEvent`.
