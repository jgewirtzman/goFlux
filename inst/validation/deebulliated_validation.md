# Validation of `goAquaFlux(diffusion.window = "deebulliated")`

Script: `deebulliated_validation.R` (run from the package root on branch
`feat/aqua-diffusive-deebulliated`; ~3 min on 6 cores). Tables:
`summary_synthetic_*.csv`, `summary_sensitivity.csv`, `results_blueflux.csv`,
`events_blueflux.csv`; per-trace synthetic results in `results_synthetic.csv.gz`.

## A. Synthetic traces (2400, seed 20260924)

Random factorial over slope {-0.5, 0, 0.05, 0.5, 2, 10} ppb/s, sigma {0.5, 2, 5} ppb
(Gaussian, 35 % AR(1) phi = 0.45), length {120, 300, 600} s, dt {1, 3, 5} s, bubbles
{0, 1, 2, 4} with steps log-uniform 5-2000 ppb, overshoot {0, 0.2, 0.5} x step with
tau {5, 12, 30} s, ramp {1, 5, 15} s, random timing (25 % first bubble < 30 s, 30 %
of multi-bubble traces with a pair 8-19 s apart), 25 % HM-curved diffusive part
(k = 0.002 s-1). Bias = (estimate - true)/|true| (traces with slope != 0), medians.

| subset | n | hit rate | FP / trace | step bias | pre NA | naive LM | pre_bubble | deebulliated |
|---|---|---|---|---|---|---|---|---|
| all traces | 2119 | 0.52 | 0.006 | 0.2 % | 32 % | 36 % | 1.6 % | 1.7 % |
| 0 bubbles | 512 | - | 0.004 | - | 0 % | 0 % | 0.3 % | 0.3 % |
| 1 bubble | 673 | 0.71 | 0.003 | -0.1 % | 37 % | 34 % | 1.7 % | 1.6 % |
| 2 bubbles | 522 | 0.53 | 0.011 | 3.6 % | 44 % | 158 % | 1.9 % | 3.8 % |
| 4 bubbles | 412 | 0.44 | 0.007 | 2.7 % | 49 % | 324 % | 11.7 % | 25.7 % |
| first bubble < 30 s | 621 | 0.46 | 0.006 | 0.2 % | 64 % | 102 % | 40 % | 7.5 % |
| first bubble > 120 s | 476 | 0.61 | 0.004 | 0.3 % | 6 % | 86 % | 1.4 % | 1.6 % |
| length 120 s | 377 | 0.37 | 0.008 | 0.8 % | 51 % | 407 % | 59 % | 17 % |
| smallest step < 3 sigma | 206 | 0.31 | 0.005 | 0.9 % | 26 % | 53 % | 7.8 % | 7.4 % |
| smallest step > 100 sigma | 400 | 0.77 | 0.007 | 0.0 % | 51 % | 289 % | 1.3 % | 1.5 % |
| pair closer than 20 s | 528 | 0.39 | 0.008 | 32 % | 46 % | 278 % | 5.1 % | 17 % |
| AR(1) noise | 764 | 0.53 | 0.005 | 0.4 % | 30 % | 44 % | 1.5 % | 1.7 % |
| HM curvature | 514 | 0.50 | 0.006 | 1.3 % | 33 % | 43 % | 0.1 % | 0.2 % |

281 traces (120 s at 5 s = 25 samples) are refused by `find.bubbles()` (>= 30
observations) and are not in the table.

**What it gets right.** Single events well above the noise are found (94 % for
one bubble with step > 30 sigma; 96-100 % for multiple bubbles more than 120 s apart)
with an unbiased settled step (median 0.0 %, fig 1), and the false-positive rate on
bubble-free traces is 0.4 % (1.6 % for AR(1) at sigma 0.5; fig 3). When every event is
found (37 % of bubble traces, mostly the crowded designs are missed) the de-ebulliated
diffusive flux has a median bias of 1.3 % and the total flux 0.2 %. The gain over the
default is where the first bubble is early: below 30 s `pre_bubble` returns NA on 64 %
of traces and is 40 % biased on the rest, `deebulliated` returns a value on 96 % with a
7.5 % median bias (fig 2); on 120 s incubations 59 % vs 17 %. Neither AR(1) wander
(fig 4) nor HM curvature changes the picture (HM is selected and its initial slope is
right). The naive whole-trace LM is 34-324 % off as soon as one bubble is present.

**Where it fails.** (1) *Closely spaced bubbles*: with gaps under 20 s the two rises are
merged into one event (hit rate 0.39, merged step +32 %), and under 60 s the second
event is often lost inside the first one's fit window (t_p + 3 tau, up to 90 s at
tau = 30) - hit 0.42-0.60 even for large steps; a missed event leaves its step in the
"de-ebulliated" trace, so the diffusive bias grows with the number of bubbles (25.7 %
median for 4 bubbles, 48 % when 4 bubbles sit within 20 s of each other) and the
summed ebullitive flux is 16 % low when any event is missed. On these traces
`pre_bubble` (when it exists) is the safer diffusive estimate. (2) *Steps at the noise
floor*: below 3 sigma only 31 % are found; the ones found are unbiased, the missed
ones bias both windows by ~7 %. (3) *Slow ramps*: 15 s ramps drop single-event
detection from 0.96 to 0.88 (0.85 at 1 Hz) because the increments per sample are
small; the found ones are still unbiased (-0.3 %). (4) *Events in the last ~10 s* are
not fitted (no post-event observations). (5) *Short incubations*: 120 s traces with
four bubbles are undetected 42 % of the time (the adaptive threshold sees the bubbles as
the norm). (6) Whatever is missed or mis-fitted goes straight into the de-ebulliated
trace: the option is only as good as `find.bubbles()`, hence `result$deebulliated`
for plotting.

**Detector settings (`summary_sensitivity.csv`, 1683 hard traces).** `bubble.window.size = 15`
(the `find.bubbles()` default; `goAquaFlux()` passes 30) raises the hit rate on multi-bubble
traces from 0.48 to 0.55, halves the 90 % quantile of the de-ebulliated diffusive error
(419 % to 216 %) and of the step error, at the cost of false positives on bubble-free traces
rising from 0.4 % to 1.4 %. `min_gap = 5`, `k = 3` and `max_reg_window = 60` change nothing
useful. I would make `bubble.window.size = 15` the `goAquaFlux()` default (matching
`find.bubbles()`), and keep `diffusion.window = "pre_bubble"` as the default because
`deebulliated` inherits every detector miss.

## B. BlueFlux floating-chamber placements

180 of the 191 placements in `output/ebullition/all_traces.rds` (campaign segmentation;
11 have < 40 samples), floating-chamber geometry 324.3 cm2 / 4.318 L, campaign flag =
point-to-point CH4 jump >= 0.10 ppm. Nine placements carry campaign jumps (41 clustered
events); the matched sample is the 113 unflagged placements from the same analyzer-days;
`goAquaFlux()` was run with both windows on all 180.

**Detection agreement.** goFlux finds 5 events on the flagged placements, all matching a
campaign event, and misses 36 campaign "events" that are not bubbles by inspection (fig 6):
20 are the individual 6 s samples of a 39 ppb/s staircase on a Picarro trace
(`BL60_P09`, threshold scaled by dt), 12 are the samples of a 25,000 ppb runaway trace the
campaign itself excluded (`CP40_P07`; goFlux finds its one real burst), two are the return
from data gaps that the segmentation bridged with a straight line (`SRS5_P06`, and `CP40_P12`
where goFlux flags it too), two are kinks on a steep concave rise (`BL60_P03`, dt 2.8 s),
and one is a burst in the last 10 s of `CP40_P02` (no post-event samples: real miss). On the
171 unflagged placements goFlux reports 13 events (10 in the matched sample), 10 of them
below the campaign's 100 ppb per-sample threshold (6-96 ppb, i.e. below what the campaign
could see) and three above it: `SRS6_P03` (a 1087 ppb rise spread over 25 samples, none
above 100 ppb: a real bubble the fixed threshold cannot catch), `LGR2 CP40_P07` (239 ppb
with a 2.9 s overshoot, looks real) and `LGR3 CP40_P27` (a 100 ppb spike-and-return at
500 s, a disturbance: false positive, effect on the diffusive slope 0.10 -> 0.09 ppb/s).

**Diffusive flux.** The two windows differ on 9 of 171 placements where both exist (all
with a goFlux event); on the rest they are identical (the pre-bubble window is the full
trace). `pre_bubble` is NA on 9 placements (bubble before 30 samples), `deebulliated` on
none. Where the campaign's own step-corrected slope exists for a matched event the
de-ebulliated value tracks it (`CP40_P06` 4.24 vs 4.81, `CP40_P10` 20.2 vs 17.1,
`SRS6_P01` 20.0 vs 14.2 ppb/s; corr 0.56 vs 0.13 for `pre_bubble`, which is NA or the
full trace on these). The disagreements are the artefact traces: `CP40_P07` (168 vs 56)
and `CP40_P12` (1.11 vs -0.08, a gap bridged by a line).

**Bottom line.** On real floating-chamber data the de-ebulliated option gives a diffusive
flux where the default gives NA and agrees with the campaign's hand-checked
step-corrected slopes; the detector is more specific than a fixed per-sample threshold
(staircases, gaps, steep ramps) and more sensitive to slow multi-sample rises, and its
weaknesses are events in the last seconds, closely spaced or crowded events, and
disturbances that look like a step. Recommend `bubble.window.size = 15` as the
`goAquaFlux()` default and plotting `result$deebulliated` before using the option.

Figures: `fig1_step_sigma.png` (detection and step bias vs step/sigma),
`fig2_timing.png` (diffusive bias vs first-bubble time; NA fraction),
`fig3_false_positives.png` (bubble-free traces), `fig4_ar1_curvature.png`,
`fig5_blueflux.png` (pre vs deebulliated slope; event counts), `fig6_worst_traces.png`
(all flagged placements and the largest unflagged goFlux events; red dotted = campaign
jumps, orange = goFlux event windows and peaks, blue = de-ebulliated trace).
