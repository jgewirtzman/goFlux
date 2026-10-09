# Precision equivalence (2026-10-08)

`empirical.prec(method = "hadamard")` replaces two copies of the per-closure
second-difference precision. `build_reference.R` computed, for the 13 bundled
closures (10 `example.closures()`, 2 `manID.G2201i`, 1 `manID.UGGA`) and two
gases (26 closure x gas values; 1 s and 5 s logging, one closure with gaps):

| implementation | source | result vs `empirical.prec()` |
|---|---|---|
| `precision_hadamard()` | fluxqc 0.2.4, commit ecc7710, `R/precision.R` | identical (precision, interval, number of second differences, lag-1 autocorrelation; tolerance 1e-12) |
| `closure_noise()` | whole_tree_flux, `data processing/goFlux_reprocessing/precision_helpers.R` | identical (same quantities); its `t_sec` equals `closure.time()` |

`tests/testthat/test-precision.equivalence.R` checks the fork against the
frozen `reference.csv` on every test run. Lag-1 autocorrelation of the second
differences (white noise: -2/3): -0.74 to -0.56 on the example closures except
`ebullition_floating` CH4 (+0.25, the bubble, by design); the G2201i (Picarro)
CO2 records give -0.34 and -0.38, i.e. the check flags them as red or
quantized noise (their precision is 0.004 ppm).

Not compared here: `ch4-data-filtering/scripts/00_setup.R` (`sigma_closure()`,
no gap handling; differs by up to ~4 % on closures with gaps) and
`analysis/precision_bakeoff/estimators.R` (comparison code, kept as is).

Duration: goFlux's native `closure.time()` was checked against whole_tree's
`goFlux_at_interval()` adapter on goFlux 0.4.0 (manID.UGGA thinned to 1, 5 and
10 s, CH4 and CO2): MDF, LM and HM fluxes, kappa and k.max identical.
