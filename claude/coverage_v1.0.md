# pipnick test coverage assessment, v1.0

**Based on commit `bdaf1450f6756abf88ad330b6b92f1aafa4b6f89`** (`bdaf145`, "Merge pull
request #15 from kbwestfall/detection"), assessed on branch `coverage`, 2026-10-29.

## Summary

The test suite (223 tests, 1 skipped) covers **85%** of pipnick's statements and
branches. Most packages are at 87–95%. The gaps that matter are not spread evenly:
they cluster in options that are documented and offered to users but never run by a
test, and in the command-line scripts.

- **Core gaps (to address now): 11 groups, listed below by risk.** The highest-risk
  gaps are silent: if the untested code were wrong, the step would finish "OK" with
  a wrong result. Examples are the sign of the mosaic `--scale calib` offsets, a
  calibration with fixed terms (`--fix`), outlier rejection in the calibration fit,
  and the adapter to astrometry.net (its 1-based pixel convention).
- **Corner cases (lower risk, for later):** error messages for invalid input,
  warnings, logging internals, maintenance tools that need the network, and the
  branches of the QA figures and the report.
- **Not a test gap:** about 110 of the 231 missed statements in `photometry/` are
  dead code (about 80 in `moffat.py`, mostly unused derivative functions, and about
  30 in `starfind.py`). It should be removed, not tested (see below).

## Method

```
COVERAGE_FILE=<scratch>/.coverage python -m pytest pipnick/tests \
    --cov=pipnick --cov-branch --cov-config=pyproject.toml --cov-report=json:<scratch>/cov.json
```

- Python 3.13.14, photutils 3.1.1, astropy 8.0.2, numpy 2.5.4, coverage 7.16.2.
- Statement **and branch** coverage, with the `[tool.coverage]` configuration in
  `pyproject.toml` (the tests and `pkg/version.py` are excluded).
- Every missed line and branch was mapped to its function, and the code was read to
  judge what the gap means. A missed line in an error message is not a missed
  algorithm.

**Classification.**

- **Core:** code that produces a product, or decides what a product contains, on a
  path users take: the default path, the options documented on the command line,
  and the degraded inputs that real nights produce (no flats in a filter, a sparse
  field). A defect would give a wrong product, often silently.
- **Corner case:** the rejection of invalid input, warnings and notes, diagnostics,
  infrastructure that only formats output, and tools used once to build packaged
  data. A defect would give a wrong message or a crash with a clear error, not a
  wrong product.

## Coverage by package

| Package | Statements | Missed | Branches | Missed | Coverage |
|---|---|---|---|---|---|
| astrometry | 680 | 48 | 202 | 20 | 92% |
| cameras | 324 | 18 | 110 | 15 | 92% |
| datamodel | 649 | 38 | 214 | 35 | 92% |
| image | 102 | 3 | 30 | 3 | 95% |
| metadata | 518 | 27 | 206 | 21 | 93% |
| mosaic | 313 | 37 | 94 | 34 | 83% |
| photometry | 1174 | 231 | 314 | 106 | 77% |
| pipelines | 1382 | 87 | 522 | 91 | 91% |
| pkg | 519 | 45 | 156 | 37 | 88% |
| qa | 228 | 14 | 56 | 22 | 87% |
| reduction | 186 | 16 | 64 | 15 | 88% |
| report | 392 | 37 | 162 | 37 | 87% |
| scripts | 381 | 293 | 34 | 23 | 24% |
| utils | 163 | 8 | 56 | 8 | 93% |
| **Total** | **7017** | **902** | **2220** | **467** | **85%** |

The lowest modules are `scripts/*` (0%, except `scriptbase` and `rgb`),
`photometry/landolt.py` (45%), `photometry/starfind.py` (60%), `pkg/logger.py`
(61%), `photometry/moffat.py` (66%), and `mosaic/prepare.py` (73%). The 0% of the
scripts is a core gap; most of the others are corner cases or dead code.

## Core gaps

Each gap gives what is untested, where, what could go wrong, and a test that would
cover it. All of these can use the synthetic data the tests already build
(`write_dithers`, `write_night`, `write_field`, the reduction's `night` fixture, and
the astrometry's `FakeSolver`), with no network.

### High risk: a defect would give a wrong product silently

**1. Mosaic flux scaling from the photometric calibration (`--scale calib`).**
`mosaic/prepare.py:149-157` is never run. The offsets are
`zp[reference] - zp[image]`; a reversed sign would scale every image by the inverse
of its correct factor, and the mosaic would still be built and reported OK.
*Test:* mosaics of images with known, different `PNCALZP` (written into the headers
of `write_dithers` images): the `INPUTS` offsets equal the differences in zero
point, with the right sign, and the fluxes of stars agree between the images after
scaling. Also test that a missing `PNCALZP` raises the documented error.

**2. Calibration with fixed terms (`--fix`).** `nightfit._resolve_fixed`
(`photometry/nightfit.py:41-55`) never applies a fixed value: not a number, not a
value per band, not `airmass=nominal`. No test calls `calibrate_all(fix=...)`.
`test_parse_fix` only parses the option strings. The fit subtracts the fixed terms
before fitting the rest (`photometry/calibrate.fit_calibration`). A mistake there,
or in matching a value to its band, would shift the zero point of every band fit
with the option. This is the usual recourse for a night with a small airmass range.
*Test:* `write_night` with known coefficients; fix `airmass` at its true value (as a
number, per band, and as `nominal`) and check that the fitted zero point and colour
term match the truth, and that the fixed term is recorded (`source` `fixed`) in the
calibration file and the headers (`PNCALEXT`).

**3. Outlier rejection in the calibration fit.** The rejection loop in
`calibrate_band` (`photometry/nightfit.py:331-332`) never rejects a measurement in
any test. The only outlier test uses `rej=100` to keep the outlier, in order to test
the leave-one-out residuals. A mistake in the loop (e.g., in the mask that is
updated) would bias the zero point whenever a standard is contaminated. That's the
case rejection exists for.
*Test:* `write_night(outlier=...)` with the default `rej`: the outlier is
`rejected` in the residual table, the coefficients are unbiased, and `n_rejected`
in the calibration file and summary is 1 per band.

**4. The astrometry.net adapter (`AstrometryNetSolver`).**
`astrometry/services.py:96-98, 132-133, 152-159, 175-182` never run, because the
tests replace the whole solver with `FakeSolver`. Untested:
- the conversion of pipnick's 0-based positions to astrometry.net's 1-based ones
  (`x + 1`, `y + 1`): a mistake would shift every astrometry.net solution by one
  pixel, then propagate through bootstrapping to every image of the setup;
- the conversion of a result into `solved` or `failed`;
- the conversion of a `TimeoutError` into `pending`, with the submission ID, for
  submission, upload and monitoring.

*Test:* an `AstrometryNetSolver` whose astroquery client (`_client`) is replaced by a
small fake that records its arguments and returns a header, an empty result, or
raises `TimeoutError(msg, submission_id)`. This needs no network, and keeps the
tests' blocking of the real service.

**5. Accepting astrometric solutions.** Two paths that decide what is written to the
headers are untested:
- *Checking an astrometry.net solution* (`astrometry/solve._check_submission`,
  `:597-601`, and its use in `handle_result`, `:639`). A solution with the wrong
  pixel scale, or too far from the pointing, should be rejected; no test returns
  one, so a bad solution could be accepted unnoticed.
- *Bootstrapping from 4 or 5 matched sources* (`solve._fit_matched:449-451`,
  `fit.wcs_from_transform:259-271`). With fewer than 6 matches, the solution is the
  reference's WCS transformed by the matched rotation and shift. That's the normal
  route for sparse fields (the default `--min_matches` is 4), but every test matches
  at least 6.

*Tests:* a `FakeSolver` that returns a solution with a wrong scale, then one far from
the pointing: the image stays unsolved, with the reason. And a sparse
`write_field`, or `max_sources` lowered, so that an image is bootstrapped from 4–5
matches: its solution agrees with the truth to within the expected error.

### Medium risk: a defect would break a common path, or degrade a product

**6. The command-line scripts.** `scripts/*.py` are at 0% (except `scriptbase` and
`rgb`), so neither their parsers (`get_parser`) nor their `main` are run. Each
`main` passes some 5–15 options to its pipeline. A misspelled keyword, a dropped
option, or a parser default that disagrees with the pipeline would only fail on a
user's command line. `pipnick_metadata`'s `main` also holds about 37 lines of
validation logic (`--split`, `--append`, `--overwrite`, existing files).
*Test:* one smoke test per script, calling `Script.main(Script.parse_args([...]))` on
the synthetic data the pipeline tests already build, as `test_rgb.py` does.
Check that it runs, that its exit status matches the verdict, and that the
non-default options reach the pipeline (e.g., `--no_qa` writes no figures).

**7. The default reduction directory.** When no `-r` is given,
`DataSet.from_file` (`datamodel/dataset.py:177-180`) sets the reduction directory
to the parent of the raw directory, and **writes it into the user's metadata
table**. `reduction_dir` (`:577-583`) has the same fallback for the scripts' logs.
That's every user's first run, but every test passes `rdxdir`.
*Test:* reduce the `night` fixture without `rdxdir`: the products are in the raw
directory's parent, the table records `rdxdir`, and its rows (including commented
ones) are unchanged.

**8. Removing an out-of-date calibration from a header.**
`photometry/apply.apply_to_headers:85-88` deletes the `PNCAL*` and `PHOTZP` keywords
from an image whose band is no longer calibrated (e.g., after a refit in which the
band failed). If this failed, the image would keep a zero point that no longer
applies, and downstream steps (`--scale calib`, the report) would use it.
*Test:* calibrate and apply; then refit with a setting that fails one band
(e.g., `min_stars` too high for it) and apply again: its images lose the keywords,
and the others keep theirs.

**9. Mosaic options and fallbacks.** These mosaic options, all offered by
`pipnick_mosaic`, are untested:
- an explicit `--reference` (`mosaic/prepare._reference:109-111`);
- the scaling by PSF fluxes, used with `--psf_match none --scale sources`
  (`prepare._source_mags:91-93`);
- the fallback from a Moffat to a Gaussian kernel when the Moffat kernel fails
  (`mosaic/psfmatch.py:113, 120, 199`; only the Gaussian kernel requested directly
  is tested);
- the rejection of outlying magnitude differences, and the inflation of the errors
  by the reduced chi-square, in the scale solution (`mosaic/scaling.py:102, 115`);
- the `interp` and `exact` reprojections (`mosaic/coadd.reproject_image:154-161`);
- the `PHOTZP`, `PNWARN` and `PNWARN2` keywords (`pipelines/mosaic._product:502-507`).

*Tests:* extend `test_mosaic_options` and the scaling tests: a named reference
gives it offset 0; PSF-flux scaling recovers the injected transparencies; a kernel
forced to fail (e.g., a target narrower than the PSF) falls back, with
`kernel` recorded as `gaussian (fallback)`; one contaminated magnitude difference
is rejected; each reprojection conserves the flux of a star.

**10. The reduction's degraded inputs.** These paths occur in real nights and are
untested:
- no bias (`pipelines/reduction.reduce_all:181`, and the science frames reduced
  without one);
- a filter without flats (`:204`);
- products left over from an earlier reduction (`:217`);
- refusing a stale source table (`ReducedImage.photometry`,
  `datamodel/image.py:280-293`; the check itself, `photometry_state`, is tested).

*Tests:* the `night` fixture with its biases (or the flats of one filter) commented
out: the verdict needs attention, the products record `PNBIAS`/`PNFLAT` `none`, and
`pipnick_status` flags `no bias`/`no flat`. A source table whose image has changed:
the calibration refuses it, with the documented error.

**11. The default radius of `--assign_mosaics`.** Without `--mosaic_radius`, the
pointings are grouped within half the field of view, computed from the camera
(`metadata/groups.assign_mosaics:224-231`). Every test passes a radius.
*Test:* `write_dithers` without `mosaic_radius`: the dithers of one field share a
label, and a field more than a field of view away does not.

## Corner cases (lower risk)

Not to address now. Grouped by kind, with the main locations.

- **Rejection of invalid input:** unknown options and malformed values
  (`pkg/settings.validate`, `pipelines/calibrate:274, 278`,
  `photometry/calibrate.check_terms`, `nightfit._resolve_fixed:41, 48-49`,
  `starfind.extract_stamps:323, 333`), missing directories (`pipelines/*: "does
  not exist"`), a reduction directory from another setup
  (`dataset.reduction_state:396, 398`).
- **Warnings and notes that do not change a product:** a band without a nominal
  extinction (`nightfit:295-296`), a small colour span under `--force`
  (`:312-313`), a negative fitted extinction (`:348`), the colour of a band not
  including it (`pipelines/calibrate:295`), frames without coordinates
  (`reduction._check_pointing:472-480`), a bias older than the measurement keywords
  (`_check_measurements:552`).
- **Metadata edge paths:** a changed reduction directory or classification on
  `--append`, no new files to append, a default table name
  (`metadata/build.py:418-499`), and frames without coordinates in
  `assign_mosaics` (`groups:247, 252`).
- **Camera checks not used by the packaged configurations:** a bad-pixel mask file
  (`cameras.unbinned_bpm:690-699`; no camera configuration sets one), and the
  temperature, empty-image and unknown-readout rejections in `assess`
  (`:827-847`). The unknown-readout rejection is worth testing once the C2 fast
  readout is confirmed, since it decides whether those frames are reduced.
- **Infrastructure:** the logger's exception hook, the recovery of the frame of a
  `warnings.warn`, and the closing of the log file (`pkg/logger.py:225-335`); the
  status and report branches for absent summaries (`report/status.py`,
  `report/records.py`, `pipelines/report.py`).
- **QA figures and displays:** the elliptical-PSF branch of `qa/psf.plot_psf`
  (`:49-54`; `test_elliptical_psf` now runs with `qa=False`), the interactive
  window of `pipnick_inspect --display` (`qa/inspect.py:66-67`, which needs a
  display), and minor branches in the other figures. A failed figure is
  recorded as needing attention, and doesn't affect a product.
- **Maintenance tools:** `photometry/landolt.build_catalog`, `_ingest`, and
  `derive_magnitudes` (`:91-221`), which build the packaged Landolt catalogue from
  VizieR, once, with the network.
- **Utilities:** degenerate polygons in `utils/masking.py`, and stretch options in
  `utils/rgb.normalize_channels`.

## Dead code

Code that no path calls is reported as missed coverage, but should be removed
rather than tested (as in Step 1 of the review). Found while reading the gaps:

- `photometry/moffat.py`, about 80 missed statements:
  - `FitMoffat2D.deriv_resid` and `FitEllipticalMoffat2D.deriv_resid`, and
    `_eval_moffat_deriv`. The Jacobian of the fit is commented out
    (`#, jac=self.deriv_resid`, `:262`).
  - `EllipticalMoffat2DModel.fit_deriv`, `input_units`, and
    `_parameter_units_for_data_units`, and `FitEllipticalMoffat2D.get_nice_phi`.
  - The `image_moments=False` branches of `guess_par`.
- `photometry/starfind.py`, about 30 missed statements: `plot_stamps`, reachable
  only through the `pdf_file` option of `build_moffat_psf_model`, which no caller
  passes; and the path and `CCDData` inputs of `as_masked_array` (the callers pass
  masked arrays).

## Notes on the tooling

- **Coverage is not run in CI.** `.github/workflows/ci_tests.yml` runs the `test`
  environments only; the tox `cov` environments exist but are not in the matrix.
  Adding one (e.g., `3.13-test-cov`) would track coverage over time.
- **An exclusion pattern never matches.** `[tool.coverage.report] exclude_lines` in
  `pyproject.toml` contains `"'def main(.*):'"`. The single quotes inside the string
  are part of the regular expression, so it never matches `def main(cls, args):`.
  This is fortunate: the scripts' `main` methods are measured, and should be. The
  pattern should be removed rather than fixed.
- **Branch coverage.** The `cov` environment measures statements only;
  `--cov-branch` (or `branch = true` in `[tool.coverage.run]`) also finds the
  untaken sides of conditions. Several gaps above (e.g., rejection that never
  rejects) show up only as branches.

## Recommended order

1. Items 1–5 (high risk), each a focused test on the existing synthetic data.
2. Item 6 (script smoke tests), which is cheap and covers the user interface.
3. Items 7–11.
4. Separately: remove the dead code, add a coverage environment to CI, and fix the
   coverage configuration.

## Follow-up: dead code (2026-10-29, branch `coverage`)

At the author's request, the dead code was dealt with before any tests were added:

- **`photometry/starfind.py`: removed.** `plot_stamps`; the `pdf_file` option of
  `psf_photometry.build_moffat_psf_model`, which called it; the path and `CCDData`
  inputs of `as_masked_array` (now an array, masked or not); and the imports
  only they used (`matplotlib`, `CCDData`, `read_reduced`).
- **`photometry/moffat.py`: commented out**, each block preceded by a note that
  the code is viable but not currently used by the pipnick pipelines: both
  `deriv_resid` methods, `FitEllipticalMoffat2D._eval_moffat_deriv`, and
  `EllipticalMoffat2DModel.fit_deriv` (the class falls back to the base class's
  `fit_deriv = None`, so an astropy fitter would estimate the derivatives
  numerically).
- **The rest of the unused code in `moffat.py`: commented out the same way.**
  `FitEllipticalMoffat2D.get_nice_phi`, and the unit hooks of
  `EllipticalMoffat2DModel` (`input_units`, `_parameter_units_for_data_units`;
  the note adds that they need `UnitsError`, whose import was removed).  In both
  `guess_par` methods, the image-centre guess of the centre and amplitude is
  commented out, and the `image_moments` argument that selected it is removed,
  so the image moments are always used (as every caller already did).
- The tests pass (223 passed, 1 skipped), with the line-length check, and the
  docs build without warnings.

## Follow-up: tests of the core gaps (2026-10-29, branch `coverage`)

Tests were added for items 1–11, each run on its own, then the full suite.

| Item | Tests |
|---|---|
| 1 | `test_mosaic.py::test_mosaic_calib_scale` |
| 2 | `test_calibrate.py::test_fixed_terms` |
| 3 | `test_calibrate.py::test_outlier_rejection` |
| 4 | `test_astrometry.py::test_astrometry_net_adapter` (with `FakeClient`) |
| 5 | `test_astrometry.py::test_check_submission` (with `DistortingSolver`), `test_fit_from_few_matches` |
| 6 | `test_scripts.py` (every script's `main`, on synthetic data) |
| 7 | `test_reduction.py::test_default_reduction_dir` |
| 8 | `test_calibrate.py::test_remove_header_calibration` |
| 9 | `test_mosaic.py::test_mosaic_reference_and_psf_fluxes`, `test_mosaic_mixed_filters`, `test_kernel_fallback`, `test_solve_offsets_rejection`, `test_reproject_methods` (and `PHOTZP` in `test_mosaic_calib_scale`) |
| 10 | `test_reduction.py::test_reduce_without_calibrations`, `test_reduce_stale_outputs`; `test_photometry.py::test_stale_photometry_guard`; `test_calibrate.py::test_stale_photometry` |
| 11 | `test_mosaic.py::test_default_mosaic_radius` |

- **Result:** the full suite passes (247 passed, 1 skipped), and the docs build
  without warnings.  Every line named for items 1–11 above is now covered;
  the total coverage is 92% (from 85%, partly from the removal of the dead
  code), and the scripts' statements are 97% covered (from 23%).
- **Defect found:** the comment of the mosaic `PHOTZP` keyword was truncated in
  every mosaic with a calibrated reference image (49 characters do not fit
  after the 20-character value field).  It is now `Reference image zero point
  (e-/s), color 0`, and `test_mosaic_calib_scale` checks it.
- **No other defects** were found.  Two test failures during the work were in
  the tests: astroquery raises its own `TimeoutError`
  (`astroquery.exceptions.TimeoutError`), which the adapter correctly catches;
  and the default mosaic radius scales with the frame size, so its test uses
  full-size frames.
