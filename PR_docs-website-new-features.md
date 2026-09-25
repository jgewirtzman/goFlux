# docs: website pages for the new features

**Branch:** `docs/website-new-features` (off master; the pages describe the functions of the
feature branches, so merge it after them)
**Files:** `website/precision.qmd` (new), `website/qualityflags.qmd` (new), `website/_quarto.yml`
(two sidebar entries), `website/manualID.qmd` (`click.peak2(gases = )`), `website/goAquaFlux.qmd`
(de-ebulliated window, one validation figure), `website/import.qmd` (`import.LI7810(dates = )`),
`website/other.qmd` (`crop.meas` aux forms, order with `goAquaFlux()`),
`website/images/deebulliated_timing.png`

The site is the Quarto project in `website/`, rendered by `.github/workflows/docs.yaml` to the
`gh-pages` branch on every push to master (the separate Qepanna/goFlux-webpage repository
publishes to quarto.pub and no longer has the goAquaFlux page, so it is not the live source).
Usage/argument sections of existing pages come from `scripts/autodoc.R` at render time; the two
new pages use plain `eval: false` code chunks on `manID.UGGA` and do not call `autodoc()`, so
they render before the feature branches are merged. There is no changelog page on the site;
the `dates =` and `crop.meas` notes went to the import and other-functions pages.
