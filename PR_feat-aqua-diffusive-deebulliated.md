# feat: de-ebulliated diffusive window in goAquaFlux()

**Branch:** `feat/aqua-diffusive-deebulliated`
**Files:** `R/goAquaFlux.diffusive.R` (new arg, `.deebulliate()` helper), `R/goAquaFlux.R` (pass-through, `diffusive_window` column), `man/goAquaFlux.Rd`, `man/goAquaFlux.diffusive.Rd`, `tests/testthat/test-goAquaFlux.deebulliated.R` (new; adds `testthat` to Suggests)

`goAquaFlux.diffusive()` fits the diffusive flux of the bubble gas on the
observations before the first bubble, so an incubation that bubbles in the
first `diffusion.minimum_window` observations gets `flux_diffusive = NA` even
though `find.bubbles()` has already fitted a model to every event. This adds
`diffusion.window = c("pre_bubble", "deebulliated")` to `goAquaFlux()` (and
`diffusive.window` to `goAquaFlux.diffusive()`). With `"deebulliated"` each
event's fitted model is subtracted from the trace (for t >= `t.step`:
`magnitude + overshoot * exp(-(t - t.peak)/tau)`, or the plain
`magnitude.step` when no re-equilibration term was retained), the samples of
the physical rise itself (from the event's detected `start` to `t.peak`, and
through `t.peak + tau` when the re-equilibration term was retained) are set
to `flag = 0`, and LM/HM are fitted on the remaining observations of the
whole incubation. The window used is returned in a new `diffusive_window`
column of `flux_summary` (`"full"`, `"pre_bubble"`, `"deebulliated"`), and
the de-ebulliated traces (with the excluded samples at `flag = 0`) come back
as a fourth list element, `result$deebulliated`, so they can be plotted.
Default behaviour is unchanged; the ebullitive flux is not touched by this
option; for gases other than the bubble gas the existing
pre-bubble/abrupt-change rule is kept because the event models were fitted on
the bubble gas.

Results: synthetic 360 s trace, slope 2 ppb/s, two bubbles (steps 400 and 250
ppb with 150/80 ppb overshoots, tau 12 s) at 60 and 200 s: de-ebulliated
diffusive flux 8.30 vs true 8.23 (0.8 %, n = 307 after excluding the two
rises) against 8.56 for the 45-obs pre-bubble window (4 %); ebullition 7.50
vs true 7.43 (0.9 %); both events recovered (400.5/150.7/11.8 and
253.3/79.9/10.8); the retained de-ebulliated samples fit a straight line of
slope 2.00. With the first bubble moved to 20 s the pre-bubble window returns
NA and the de-ebulliated fit still gives 8.3. On a real floating-chamber burst
closure (step 1089 ppb at 23 s, 278 obs) the default gives NA (24 pre-bubble
obs) and the de-ebulliated fit gives 22.6 (HM, n = 251 after flagging out the
27 samples from the detected start at 23 s to the peak at 49.8 s; 21.1 when
those samples were left in; the LM slope matches the raw post-bubble slope),
with the ebullitive flux unchanged at 21.08. A fuller validation (synthetic
grid and BlueFlux floating-chamber closures) is in `inst/validation/`.
