# Stage Taxonomy and Type Hierarchy

## Overview

Every `MatrixProjectionModel` carries a vector of `StageClass` annotations
classifying each stage as **active**, **propagule**, or **dormant**. This
vignette walks through the stage classification API, the constructors that
consume / produce stage metadata, and the abstract type hierarchy
classifying models by structure (Leslie vs. Lefkovitch), density
dependence, and stochasticity.

## Setup

```@example mpm
using MatrixProjectionModels
using LinearAlgebra
```

## A reference model

We use a 4-stage model with reproduction concentrated in adults:

```@example mpm
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

## The `StageClassType` enum

Stage roles are encoded by an enum with three members.

```@example mpm
println("StageClassType values: ", instances(StageClassType))
println("ActiveStage    = ", ActiveStage)
println("PropaguleStage = ", PropaguleStage)
println("DormantStage   = ", DormantStage)
```

`StageClass` wraps the enum together with an author label and a within-author
index.

```@example mpm
sc_active = StageClass(ActiveStage,    "Caswell", 1)
sc_prop   = StageClass(PropaguleStage, "Caswell", 2)
sc_dorm   = StageClass(DormantStage,   "Caswell", 3)
println("typeof(sc_active).name.name = ", typeof(sc_active).name.name)
```

## Automatic classification

`standard_stages` infers the classification from the structure of the model.
Stages that **receive** reproduction but do not themselves reproduce are
typically tagged `PropaguleStage`; otherwise they are `ActiveStage`. Dormant
stages must be supplied explicitly.

```@example mpm
classes = standard_stages(mpm.F)
for (i, c) in enumerate(classes)
    println("  stage $i (", mpm.stage_names[i], ") → ", c)
end
```

`repro_stages` returns the BitVector of reproductive columns:

```@example mpm
println("reproductive stages = ", repro_stages(mpm.F))
```

## Renaming stages

`name_stages` returns a renamed copy:

```@example mpm
mpm2 = name_stages(mpm, ["juv", "sm_ad", "lg_ad", "seedbank"])
println("renamed stage_names = ", mpm2.stage_names)
```

## Mature stage distribution

`mature_distrib` returns the fraction of individuals that reach each
reproductive stage starting from each non-reproductive stage.

```@example mpm
md = mature_distrib(mpm.U; repro_stages = repro_stages(mpm.F))
println("mature_distrib (per starting stage) = ", round.(md, digits=4))
```

## Iterating the model

```@example mpm
sol = solve(MPMProblem(mpm, [10.0, 5.0, 1.0, 1.0], (0, 50)))
println("typeof(sol).name.name              = ", typeof(sol).name.name)
println("sol isa MPMSolution                = ", sol isa MPMSolution)
println("sol isa AbstractProjectionSolution = ", sol isa AbstractProjectionSolution)
```

## Structure types

`AbstractMPMStructure` is the dispatch tag set for matrix structure. Two
canonical members are `LeslieMPM` (age-based) and `LefkovitchMPM`
(stage-based).

```@example mpm
println("LeslieMPM     <: AbstractMPMStructure              = ",
        LeslieMPM <: AbstractMPMStructure)
println("LefkovitchMPM <: AbstractMPMStructure              = ",
        LefkovitchMPM <: AbstractMPMStructure)
println("AbstractMPMStructure <: AbstractProjectionStructure = ",
        AbstractMPMStructure <: AbstractProjectionStructure)
```

Density and stochasticity tags share parallel hierarchies.

```@example mpm
println("DensityIndependent() isa AbstractDensityDependence = ",
        DensityIndependent() isa AbstractDensityDependence)
```

The `AbstractStochasticity` supertype covers stochastic dispatch tags
(deterministic / stochastic backends).

```@example mpm
println("AbstractStochasticity is an abstract type: ",
        isabstracttype(AbstractStochasticity))
```

## Summary

- `StageClassType` is the enum vocabulary for stage roles, with members
  `ActiveStage`, `PropaguleStage`, `DormantStage`.
- `standard_stages` / `repro_stages` infer the classification; `name_stages`
  is the rename helper; `mature_distrib` projects the asymptotic
  contribution onto reproductive stages.
- `MPMSolution`, `AbstractMPMStructure`, `LeslieMPM`, `LefkovitchMPM`,
  `AbstractDensityDependence`, `DensityIndependent`,
  `AbstractStochasticity`, and `AbstractProjectionStructure` /
  `AbstractProjectionSolution` are the core dispatch types.
