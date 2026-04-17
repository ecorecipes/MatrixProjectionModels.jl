# Stage Taxonomy and Type Hierarchy
Simon Frost

## Overview

Every `MatrixProjectionModel` carries a vector of `StageClass`
annotations classifying each stage as **active**, **propagule**, or
**dormant**. This vignette walks through the stage classification API,
the constructors that consume / produce stage metadata, and the abstract
type hierarchy classifying models by structure (Leslie vs. Lefkovitch),
density dependence, and stochasticity.

## Setup

``` julia
using MatrixProjectionModels
using LinearAlgebra
```

## A reference model

We use a 4-stage model with reproduction concentrated in adults:

``` julia
U = [0.0  0.0  0.0  0.0;
     0.5  0.4  0.0  0.6;
     0.0  0.3  0.5  0.0;
     0.0  0.0  0.2  0.0]
F = [0.0  0.0  3.0  0.0;
     0.0  0.0  0.0  0.0;
     0.0  0.0  0.0  0.0;
     0.0  0.0  0.5  0.0]

mpm = MatrixProjectionModel(U, F;
    stage_names = [:juvenile, :small_adult, :large_adult, :seedbank])
println("typeof(mpm).name.name = ", typeof(mpm).name.name)
```

    typeof(mpm).name.name = MatrixProjectionModel

## The `StageClassType` enum

Stage roles are encoded by an enum with three members.

``` julia
println("StageClassType values: ", instances(StageClassType))
println("ActiveStage    = ", ActiveStage)
println("PropaguleStage = ", PropaguleStage)
println("DormantStage   = ", DormantStage)
```

    StageClassType values: (ActiveStage, PropaguleStage, DormantStage)
    ActiveStage    = ActiveStage
    PropaguleStage = PropaguleStage
    DormantStage   = DormantStage

`StageClass` wraps the enum together with an author label and a
within-author index.

``` julia
sc_active = StageClass(ActiveStage,    "Caswell", 1)
sc_prop   = StageClass(PropaguleStage, "Caswell", 2)
sc_dorm   = StageClass(DormantStage,   "Caswell", 3)
println("typeof(sc_active).name.name = ", typeof(sc_active).name.name)
```

    typeof(sc_active).name.name = StageClass

## Automatic classification

`standard_stages` infers the classification from the structure of the
model. Stages that **receive** reproduction but do not themselves
reproduce are typically tagged `PropaguleStage`; otherwise they are
`ActiveStage`. Dormant stages must be supplied explicitly.

``` julia
classes = standard_stages(mpm.F)
for (i, c) in enumerate(classes)
    println("  stage $i (", mpm.stage_names[i], ") → ", c)
end
```

      stage 1 (juvenile) → PropaguleStage
      stage 2 (small_adult) → ActiveStage
      stage 3 (large_adult) → ActiveStage
      stage 4 (seedbank) → PropaguleStage

`repro_stages` returns the BitVector of reproductive columns:

``` julia
println("reproductive stages = ", repro_stages(mpm.F))
```

    reproductive stages = Bool[0, 0, 1, 0]

## Renaming stages

`name_stages` returns a renamed copy:

``` julia
mpm2 = name_stages(mpm, ["juv", "sm_ad", "lg_ad", "seedbank"])
println("renamed stage_names = ", mpm2.stage_names)
```

    renamed stage_names = [:juv, :sm_ad, :lg_ad, :seedbank]

## Mature stage distribution

`mature_distrib` returns the fraction of individuals that reach each
reproductive stage starting from each non-reproductive stage.

``` julia
md = mature_distrib(mpm.U; repro_stages = repro_stages(mpm.F))
println("mature_distrib (per starting stage) = ", round.(md, digits=4))
```

    mature_distrib (per starting stage) = [0.0, 0.0, 1.0, 0.0]

## Iterating the model

``` julia
sol = solve(MPMProblem(mpm, [10.0, 5.0, 1.0, 1.0], (0, 50)))
println("typeof(sol).name.name              = ", typeof(sol).name.name)
println("sol isa MPMSolution                = ", sol isa MPMSolution)
println("sol isa AbstractProjectionSolution = ", sol isa AbstractProjectionSolution)
```

    typeof(sol).name.name              = MPMSolution
    sol isa MPMSolution                = true
    sol isa AbstractProjectionSolution = true

## Structure types

`AbstractMPMStructure` is the dispatch tag set for matrix structure. Two
canonical members are `LeslieMPM` (age-based) and `LefkovitchMPM`
(stage-based).

``` julia
println("LeslieMPM     <: AbstractMPMStructure              = ",
        LeslieMPM <: AbstractMPMStructure)
println("LefkovitchMPM <: AbstractMPMStructure              = ",
        LefkovitchMPM <: AbstractMPMStructure)
println("AbstractMPMStructure <: AbstractProjectionStructure = ",
        AbstractMPMStructure <: AbstractProjectionStructure)
```

    LeslieMPM     <: AbstractMPMStructure              = true
    LefkovitchMPM <: AbstractMPMStructure              = true
    AbstractMPMStructure <: AbstractProjectionStructure = true

Density and stochasticity tags share parallel hierarchies.

``` julia
println("DensityIndependent() isa AbstractDensityDependence = ",
        DensityIndependent() isa AbstractDensityDependence)
```

    DensityIndependent() isa AbstractDensityDependence = true

The `AbstractStochasticity` supertype covers stochastic dispatch tags
(deterministic / stochastic backends).

``` julia
println("AbstractStochasticity is an abstract type: ",
        isabstracttype(AbstractStochasticity))
```

    AbstractStochasticity is an abstract type: true

## Summary

- `StageClassType` is the enum vocabulary for stage roles, with members
  `ActiveStage`, `PropaguleStage`, `DormantStage`.
- `standard_stages` / `repro_stages` infer the classification;
  `name_stages` is the rename helper; `mature_distrib` projects the
  asymptotic contribution onto reproductive stages.
- `MPMSolution`, `AbstractMPMStructure`, `LeslieMPM`, `LefkovitchMPM`,
  `AbstractDensityDependence`, `DensityIndependent`,
  `AbstractStochasticity`, and `AbstractProjectionStructure` /
  `AbstractProjectionSolution` are the core dispatch types.
