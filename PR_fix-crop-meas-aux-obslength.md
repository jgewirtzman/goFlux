# fix: crop.meas(max.obs.length = "aux") always failed with a NULL join

**Branch:** `fix/crop-meas-aux-obslength` (one line in `R/crop.meas.R`)

In `crop.meas()` the per-UniqueID branch for `max.obs.length` assigned the
auxfile columns to `unique_max.obs.length` but then joined
`unique_obs.length`, which is initialised to `NULL` in the "variables without
binding" line, so the documented `max.obs.length = "aux"` form errored with
"`x` and `y` must share the same src ... `y` is NULL" (the scalar form worked).

```r
data(manID.UGGA)
aux <- data.frame(UniqueID = unique(manID.UGGA$UniqueID), max.obs.length = 120)
crop.meas(manID.UGGA, auxfile = aux, max.obs.length = "aux")   # errored before, works now
```

No documentation change; behaviour of every other argument unchanged.
