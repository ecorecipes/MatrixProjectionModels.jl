# COMPADRE Integration

Interface to the [COMPADRE Plant Matrix Database](https://compadre-db.org/) and COMADRE Animal Matrix Database. These functions are provided via a package extension and require loading the appropriate backend.

## Database Type

```@docs
AbstractCompadreDB
```

## Loading & Saving

```@docs
cdb_load
cdb_save
cdb_fetch
```

## Matrix Access

```@docs
cdb_matA
cdb_matU
cdb_matF
cdb_matC
```

## Metadata & Filtering

```@docs
cdb_metadata
cdb_id
cdb_flag
cdb_subset
```

## Manipulation

```@docs
cdb_collapse
cdb_rbind
cdb_flatten
cdb_build_cdb
```

## Summary Statistics

```@docs
mpm_mean
mpm_sd
mpm_median
```

## Stage Properties

```@docs
mpm_has_prop
mpm_has_active
mpm_has_dorm
mpm_first_active
```
