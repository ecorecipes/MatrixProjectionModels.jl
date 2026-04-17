# COMPADRE/COMADRE Database Interface
Simon Frost

## Overview

The [COMPADRE](https://compadre-db.org/Data/Compadre) plant matrix
database and its animal counterpart
[COMADRE](https://compadre-db.org/Data/Comadre) contain hundreds of
curated matrix population models. This vignette documents the function
surface that `MatrixProjectionModels.jl` exposes for working with the
database via the **`MatrixProjectionModelsCompadreExt`** package
extension (loaded automatically when `CSV`, `DataFrames`, `Downloads`,
and `RData` are all present).

> **Status note.** The extension currently provides type definitions and
> a stable function surface, but RData parsing and several accessors are
> still being implemented. Calling many of these functions today will
> raise `error("... not yet fully implemented")`. The function surface
> below is the public API the extension is converging on.

## Setup

``` julia
using MatrixProjectionModels
```

## Type and method introspection

The abstract type `AbstractCompadreDB` is defined in the main module so
downstream packages can dispatch on it without forcing the extension to
load.

``` julia
println("AbstractCompadreDB is abstract: ", isabstracttype(AbstractCompadreDB))
```

    AbstractCompadreDB is abstract: true

Generic stubs for every COMPADRE function are exported from the main
package — they only acquire methods once the extension is loaded.

``` julia
api_funcs = [
    cdb_fetch, cdb_load, cdb_save,
    cdb_matA, cdb_matU, cdb_matF, cdb_matC,
    cdb_metadata, cdb_id, cdb_flag,
    cdb_collapse, cdb_rbind, cdb_flatten, cdb_subset, cdb_build_cdb,
    mpm_mean, mpm_sd, mpm_median,
    mpm_has_prop, mpm_has_active, mpm_has_dorm, mpm_first_active,
]
println("number of COMPADRE API functions exported = ", length(api_funcs))
println()
for f in api_funcs
    println("  ", rpad(string(nameof(f)), 22), " methods = ", length(methods(f)))
end
```

    number of COMPADRE API functions exported = 22

      cdb_fetch              methods = 0
      cdb_load               methods = 0
      cdb_save               methods = 0
      cdb_matA               methods = 0
      cdb_matU               methods = 0
      cdb_matF               methods = 0
      cdb_matC               methods = 0
      cdb_metadata           methods = 0
      cdb_id                 methods = 0
      cdb_flag               methods = 0
      cdb_collapse           methods = 0
      cdb_rbind              methods = 0
      cdb_flatten            methods = 0
      cdb_subset             methods = 0
      cdb_build_cdb          methods = 0
      mpm_mean               methods = 0
      mpm_sd                 methods = 0
      mpm_median             methods = 0
      mpm_has_prop           methods = 0
      mpm_has_active         methods = 0
      mpm_has_dorm           methods = 0
      mpm_first_active       methods = 0

## Loading and IO

Once the extension is loaded:

``` julia
db_compadre = cdb_fetch(:compadre; save = true, destination = "compadre.RData")
db_comadre  = cdb_fetch(:comadre)
db = cdb_load("compadre.RData")
cdb_save(db, "filtered.RData")
```

## Accessors

``` julia
matA = cdb_matA(db)            # Vector{Matrix} of A = U+F+C matrices
matU = cdb_matU(db)            # Vector{Matrix} of survival/growth blocks
matF = cdb_matF(db)            # reproduction
matC = cdb_matC(db)            # clonality
meta = cdb_metadata(db)        # the full DataFrame of metadata
ids  = cdb_id(db)              # vector of MPM IDs
flag = cdb_flag(db)            # quality flag DataFrame
```

## Operations

``` julia
db_clean = cdb_subset(db, db -> db.OrganismType .== "Tree")
db_join  = cdb_rbind(db_clean, another_db)
db_flat  = cdb_flatten(db)                         # explode list-columns
db_coll  = cdb_collapse(db, target_stages)         # collapse to fewer stages
db_built = cdb_build_cdb(metadata_df, matrices)    # construct from parts
```

## Statistics on collections of MPMs

``` julia
A_mean   = mpm_mean(matA)      # entrywise mean across replicates
A_sd     = mpm_sd(matA)        # entrywise standard deviation
A_median = mpm_median(matA)    # entrywise median
```

## Stage-class queries

``` julia
mpm_has_prop(db)      # BitVector: which entries have a propagule stage?
mpm_has_active(db)    # which have any active stage?
mpm_has_dorm(db)      # which have a dormant stage?
mpm_first_active(db)  # index of the first active stage in each entry
```

## Single-MPM transformations

Three helpers operate on a single `MatrixProjectionModel` (independently
of the database) to standardise it for cross-study comparison:

``` julia
mpm_standardize(mpm)                 # canonical orientation / sign conventions
mpm_rearrange(mpm, perm)             # permute stages
mpm_collapse(mpm, target_stages)     # collapse stages by aggregation
```

These are the building blocks that `cdb_collapse` calls per-row when
operating on a full database.

## Summary

- `AbstractCompadreDB` is the dispatch tag; concrete `CompadreDB` is
  defined inside `MatrixProjectionModelsCompadreExt`.
- IO: `cdb_fetch`, `cdb_load`, `cdb_save`.
- Accessors: `cdb_matA`/`U`/`F`/`C`, `cdb_metadata`, `cdb_id`,
  `cdb_flag`.
- Operations: `cdb_subset`, `cdb_rbind`, `cdb_flatten`, `cdb_collapse`,
  `cdb_build_cdb`.
- Aggregate stats: `mpm_mean`, `mpm_sd`, `mpm_median`.
- Stage-class queries: `mpm_has_prop`, `mpm_has_active`, `mpm_has_dorm`,
  `mpm_first_active`.

End-to-end working examples will be added when the extension’s RData
parsing layer is complete.
