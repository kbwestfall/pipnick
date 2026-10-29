# pipnick code review, v1.0

- **Branch:** `review`
- **Base commit:** `1e959c5` (main, after the `mosaic` merge)
- **Date:** 2026-10-26
- **Size:** about 16,400 lines of package code and 4,400 lines of tests.

**What was reviewed.** All of the code reachable from the eight console
scripts:

- `pipnick_version`
- `pipnick_metadata`
- `pipnick_reduce`
- `pipnick_astrometry`
- `pipnick_photometry`
- `pipnick_calibrate`
- `pipnick_mosaic`
- `pipnick_rgb`

**How it was reviewed.**

- **Unused code** was found by tracing imports and calls from the scripts.
  Everything else in this review ignores it (section 1).
- **The conventions audit** combined AST and tokenize scripts with reading of
  every flagged site.
- **The architecture and logging assessments** come from reading the
  pipelines and the modules they call.

Every finding was checked against the code. Items still in doubt are marked
*uncertain*. Line numbers refer to the base commit.

This review proposes changes only; nothing has been changed.

---

## Summary

### 1. Unused code → [details](#1-unused-code)

About 1,000 lines are not reachable from any script.

- **Unused modules:**
  - `utils/graphs.py`
  - `utils/log.py`
  - `utils/display_fits.py`, which `pipnick_reduce` imports but never calls
  - `utils/fits_class.py` (`Fits_Simple`)
- **Unused functions in `photometry/fit.py`:** three of the four functions.
- **Smaller dead items:** a few methods and helpers, plus commented-out code
  in `utils/masking.py` and `scripts/__init__.py`.
- **Unused dependencies:** `pandas`, `IPython` and `requests` are never
  imported. `loess` (which can't be redistributed) is used only by the unused
  `utils/graphs.py`.

Removing the unused code would let four dependencies go.

### 2. Conventions → [details](#2-conformance-with-agentsmd)

The code mostly conforms to `AGENTS.md`:

- no intentional NaN;
- no test asserts without messages;
- no `os.path`;
- no nested functions outside tests;
- no docstrings with a backslash missing the `r` prefix.

The violations are concentrated in four areas:

- **Wrapping style:** about 30 sites in the package and 30 in the tests,
  mostly hanging indents and trailing commas.
- **Private imports across modules (the new rule):** 3 import statements and
  7 names. The tests use 32 more.
- **Docstrings:** about a dozen missing roles or bare package names.
- **Small items:** a few, such as a truthiness check and late imports in
  `__init__.py`.

All of them are mechanical to fix.

### 3. Modularity → [details](#3-modularity)

- **Scripts:** thin and consistent.
- **Pipelines:** each mixes five concerns: finding the inputs, rerun
  bookkeeping, the science strategy, file I/O, and QA and summaries.
- **Size:** four pipeline modules and `metadata.py` each have more than
  1,000 lines.
- **Coupling:** the pipelines import from one another, and
  `pipelines/reduction.py` has become the data model of the whole package.
  `load_setup`, `plan_products` and `table_values` are used by every step.

Proposal: a `datamodel` package, an `image` package (input/output,
background, detection), a `qa` package, and slim pipeline drivers that never
import another pipeline.

### 4. Redundancy → [details](#4-redundancy-and-consolidation)

| Function | Separate implementations |
|---|---|
| Source detection | 4 (used in 5 places) |
| Background estimators | 4 |
| Source matchers | 4 |
| Hash helpers | 11 (3 byte-identical) |
| Image collectors | 4 |
| Summary writers | 4 |
| Weighted least-squares solvers | 3 |
| Rerun mechanisms | 5 |

Two of these are inconsistencies, not just duplication:

- **FWHM of an elliptical PSF:** the photometry uses the mean of the axes,
  but the mosaic uses the major axis.
- **Missing values:** some outputs use `-1` sentinels, while others mask.

**Your example is feasible.** One source-detection function can replace the
four implementations. Each step keeps running it whenever it needs to, with
its own settings, or applies its own cuts to the resulting catalogue:

- the astrometry: a high threshold and the brightest sources;
- the PSF model: well-detected, isolated stars;
- the photometry: as complete a census as possible.

One caveat: the astrometry.net cache is keyed on the detected positions. The
shared function must therefore reproduce the astrometry's current detections
exactly, which a test can check, or previously solved images could be
resubmitted.

### 5. Object orientation → [details](#5-object-orientation)

The steps pass around dictionaries (products, images, settings) and a
five-element tuple returned by `load_setup`. Proposed classes:

- `Product` and `DataSet` (a per-setup table and its reductions);
- `ReducedImage`;
- `SourceCatalog`;
- `PSFModel`, which also fixes the FWHM inconsistency;
- `Calibration` (the TOML file);
- frozen settings dataclasses;
- a `PipelineStep` base class with a common `run` loop.

The files on disk (header keywords, ECSV metadata, TOML) would not change.
The migration can go one class at a time, with the existing functions kept
as thin wrappers.

### 6. Log messages → [details](#6-log-messages)

- **Already clear:** astrometry and calibration explain their decisions, with
  reasons and suggestions.
- **Weakest:** reduction gives no progress, no closing summary and no verdict
  on its measurements (bias, read noise, flats, sky). It writes about 400
  lines a night in which most lines don't name the frame.
- **Across all steps:**
  - no pipnick version is recorded;
  - nothing is logged to a file unless `--log_file` is given;
  - no step ends with a verdict and the next command to run.

The proposals are a run header, `[n/N]` progress, checks against expected
ranges, closing verdicts and next-step hints. Section 6 gives the proposed
message text for about 25 locations.

### 7. Inspecting the results → [details](#7-evaluating-and-interrogating-the-results)

Three kinds of output exist:

- summary tables for four of the six steps;
- optional QA figures for four steps;
- full provenance in the image headers.

Missing:

- a single view of a night across the steps;
- QA figures by default, and for up-to-date or failed products;
- a summary or QA from the reduction step;
- a common status vocabulary;
- a way to trace one image or one star through the steps, or to see why a
  source was removed.

Proposals:

- `pipnick_status`, a table of the night with one row per image;
- `pipnick_report`, a static HTML report;
- `pipnick_inspect`, to examine one image or one star;
- default logs and run records;
- expected ranges for the camera, so the checks can say whether a value is
  normal.

### Suggested order of work → [details](#8-suggested-order-of-work)

1. Remove the unused code and dependencies, and fix the convention items.
   **Complete** (2026-10-26, branch `review`; see the
   [progress log](#progress-log)).
2. Consolidate the pure helpers: hashes, paths, least squares, matching.
   **Complete** (2026-10-27, branch `consolidate_helpers`; see the
   [progress log](#progress-log)).
3. Add `DataSet`, `ReducedImage` and `PSFModel`.
   **Complete** (2026-10-28, branch `common_objects`; see the
   [progress log](#progress-log)).
4. Improve the logging and add `pipnick_status`.
   **Complete** (2026-10-28, branch `status`; see the
   [progress log](#progress-log)).
5. Split the large modules, and add `SourceCatalog`, `Summary` and the `qa/`
   package.
6. Move every step to the common detection function, then add
   `PipelineStep`.

---

## 1. Unused code

The analysis followed imports and calls from the console scripts, then
checked each candidate with `grep`. Tests don't count as uses.

### 1.1 Unused modules

| Module | Lines | Status |
|---|---|---|
| `utils/graphs.py` | 267 | Never imported.  The only user of `loess` (not redistributable; see `licenses/README.rst`). |
| `utils/log.py` | 19 | Never imported (`log_astropy_table`). |
| `utils/display_fits.py` | 102 | Imported only inside `scripts/reduce.py:42`; the `--display` option and the call are commented out (`:32`, `:51-52`). |
| `utils/fits_class.py` | 369 | Used only by `display_fits.py` and `graphs.py` (`Fits_Simple`, `unzip_directories`). |

### 1.2 Unused functions and members in modules that are used

| Location | Item | Note |
|---|---|---|
| `photometry/fit.py` | `get_param_list`, `get_apertures`, `plot_sources` | Only `psf_plot` is used, by the photometry QA. |
| `photometry/moffat.py` | `FitEllipticalMoffat2D.get_orig_phi` | No callers. |
| `pkg/logger.py` | `clear_text_color` | No callers. |
| `pkg/exceptions.py` | `PipnickCodingError` | Never raised. |
| `photometry/landolt.py` | `build_catalog`, `_ingest` | Not used by the scripts, but needed to regenerate the packaged catalog.  **Keep**, as a maintenance tool. |
| `utils/masking.py` | lines 164–285 | Commented-out code. |
| `scripts/__init__.py` | — | Commented-out JSON logging configuration. |
| `scripts/reduce.py:42` | `display_many_nickel` import | Dead import. |

### 1.3 Unused dependencies

These are declared in `pyproject.toml` but never imported by package code:

- `pandas`
- `IPython`
- `requests`
- `loess`, which is imported only by the unused `utils/graphs.py`.

Removing `loess` also removes the packaging problem described in
`licenses/README.rst`. `scikit-image` is used only by
`photometry/starfind.py`.

**Recommendation:**

- Delete the modules in 1.1, the functions in 1.2 (except `build_catalog`)
  and the commented-out code.
- Drop the four dependencies.
- Update the API docs and `licenses/README.rst` to match.

Commented-out code (e.g. `scripts/reduce.py`) and the stale usage docs
(`docs/basics.rst`, `docs/example.rst`) refer to some of these items, so
clean them up together.

---

## 2. Conformance with `AGENTS.md`

The unused code from section 1 is excluded.

| # | Convention | Violations |
|---|---|---|
| 1 | No intentional NaN | 0 |
| 2 | Explicit condition checks | 1 |
| 3 | Test asserts have messages | 0 |
| 4 | Wrapping style | about 31 in package code; about 30 in tests |
| 5 | Lines ≤ 99 characters | 1 |
| 6 | No nested functions | 0 |
| 7 | Few module globals | 1 (uncertain) |
| 8 | Imports at the top | 3 (`pipnick/__init__.py`) |
| 9 | No non-standard `import ... as` | 2 (tests) |
| 10 | Docstrings with a backslash are r-strings | 0 |
| 11 | numpydoc with explicit roles | about 12 |
| 12 | TOML configuration | 0 |
| 13 | Attribution of adapted code | 1 (uncertain) |
| 14 | Sorted directory listings | 1 (uncertain, test) |
| 15 | `pathlib`, not `os.path` | 0 |
| 16 | No private names across modules | 7 names in 3 statements (plus 32 in tests) |

### 2.1 Explicit condition checks

- `pipelines/mosaic.py:626` tests a dictionary (or `None`) for truth:
  `image['phot'].meta.get('psf')`.
  - It is inconsistent with `_psf` (`:257`), which checks `fwhm_pix <= 0`.
  - Replace it with an explicit check. This is also the place where `-1`
    sentinels are written into the `INPUTS` table (see 4.6).

### 2.2 Wrapping style

**Contents begin on the opening line, then a hanging indent.** Package code:

- `pipelines/astrometry.py:904, 909`
- `pipelines/calibrate.py:393, 1220, 1583`
- `pipelines/mosaic.py:640`
- `pipelines/photometry.py:222, 533`
- `scripts/version.py:16, 30`
- `photometry/landolt.py:161, 169`

For example, `calibrate.py:393`:

```python
for name in ['sep', 'airmass', ..., 'mag_inst',
             'mag_inst_err', 'err']:
```

Tests:

- `test_astrometry.py:272, 681, 744`
- `test_calibrate.py:313, 315, 317, 324, 764, 767`
- `test_mosaic.py:84, 169, 338, 340, 505`
- `test_reduction.py:105, 151`

**Trailing comma after the last element.** Package code:

- `cameras/cameras.py:154, 886`
- `metadata.py:106` (the `metadata_columns` dictionary)
- `photometry/moffat.py:817, 831`
- `pkg/exceptions.py:16`
- `pkg/logger.py:110, 164, 369`

Tests:

- `test_cameras.py:23, 41, 171, 224, 262, 291`
- `test_metadata.py:26, 42, 58`

**One element per line where several would fit.**

- `metadata.py:88` (`metadata_columns`)
- `cameras/cameras.py:872`
- `pkg/exceptions.py:11`
- `pkg/logger.py:43, 105, 159, 365, 397`
- `photometry/moffat.py:815, 826`

A dictionary with one entry per line may be a deliberate choice for
readability. If so, record it in `AGENTS.md` as an exception.

**Backslash continuation.** `astrometry/services.py:280`.

**Uncertain.** These conditions open a call at the end of the line rather
than using `if (`:

- `pipelines/mosaic.py:913`
- `pipelines/astrometry.py:584`
- `cameras/cameras.py:703, 708`

`pkg/version.py` is generated, so it is excluded.

### 2.3 Other items

**Line length.** `photometry/moffat.py:728` is 117 characters (a URL in a
docstring). Wrap it with an RST hyperlink target.

**Globals (uncertain).** `metadata.mosaic_columns` is used by one function
in package code, but it is part of the documented interface (`__all__`,
`:data:` references) and is used by the tests. Keeping it is defensible.

**Imports.** `pipnick/__init__.py:13, 14, 20` import after
`__version__ = version`. Two of them carry `# noqa: E402`. Reorder them, or
document the exception.

**Aliases.** These test imports are effectively aliases:

- `tests/test_calibrate.py:19`: `from pipnick.pipelines import calibrate as pipeline`
- `tests/test_mosaic.py:24`: `from pipnick.pipelines import mosaic as pipeline`

Import `pipnick.pipelines.calibrate` directly, or use a qualified name.

**Docstrings.**

- `utils/combiner.py:39` uses the default role, `` `numpy.dtype` ``.
- `scripts/scriptbase.py:81` has `argparse.HelpFormatter` without a role.
- Bare package names in prose should be ``` ``literals`` ```:
  - `astrometry/fit.py:4, 25, 27`
  - `astrometry/services.py:193`
  - `cameras/cameras.py:268, 286, 723`
  - `utils/masking.py:94, 95, 111`
  - `tests/conftest.py:6`
- `pipelines/reduction.py:447, 582` cross-reference the private
  `pipnick.metadata._filter_name`. This goes away when it is made public
  (2.4).

**Attribution (uncertain).** `utils/masking.py` has no module docstring.
`polygon_winding_number` cites Numerical Recipes, and appears to match
PypeIt's `pypeit/utils.py`. If it was adapted, add the PypeIt attribution.

**Sorted listings (uncertain, test).** `tests/test_astrometry.py:874` uses
`next(...glob('*.ecsv'))`, which picks an arbitrary one of four files.

### 2.4 Private names across modules (the new rule)

| Importer | Private names | Proposed public name |
|---|---|---|
| `pipelines/calibrate.py:53` | `photometry._image_hash`, `_wcs_hash`, `_select` | `image_hash`, `wcs_hash` (move to `pkg/provenance.py`; section 4.4); `select_products` |
| `pipelines/mosaic.py:61` | `photometry._image_hash`, `_wcs_hash` | as above |
| `pipelines/reduction.py:29` | `metadata._filter_name`, `_safe_name` | `filter_name`, `safe_name` |

The tests use 32 private names, which the new rule allows. The most common
are:

- `reduction._safe_name`, `_flat_name`, `_science_name`, `_group_indices`
- `astrometry._remove_wcs`
- `calibrate._toml_value`, `_leave_one_out`
- the monkeypatched `GaiaCatalog._query`

Several of them (`_science_name`, `_flat_name`, `_group_indices`) would
become public anyway under the data model in section 5.

---

## 3. Modularity

### 3.1 What works

- The scripts are thin: `main` parses the arguments and calls one `*_all`
  function, with the imports kept inside `main`.
- The algorithm packages (`astrometry/`, `photometry/`, `mosaic/`,
  `cameras/`) are mostly free of pipeline concerns.

### 3.2 Problems

1. **The pipelines depend on each other.**
   - Every step imports `load_setup` and `plan_products` from
     `pipelines/reduction.py`:
     - `photometry.py:33`
     - `astrometry.py:63`
     - `calibrate.py:54`
     - `mosaic.py:62`
   - Photometry and mosaic also import `table_values`.
   - Calibrate and mosaic import private helpers from
     `pipelines/photometry.py`.

   `reduction.py` is in effect the package's data model.
2. **Each pipeline mixes five concerns.** For example,
   `pipelines/astrometry.py` (1,578 lines) contains:
   - image state (`AstrometryImage`, `_Session`);
   - header input/output (`_remove_wcs`, `apply_wcs`);
   - the astrometry.net cache (`_anet_*`, `_restore_solution`);
   - the solving strategy (`group_by_pointing` … `_solve_group`);
   - the submission (`_submit`, `_handle_result`);
   - reporting (`_summarize`, `_write_qa`, `_report_dry_run`).

   The calibration (1,660 lines), reduction (1,362), mosaic (1,061) and
   metadata (1,310) modules are similar.
3. **Code in the wrong place.**
   - The calibration's command-line parsers (`parse_fix`,
     `parse_assignments`, `default_colors`) are in the pipeline module but
     only used by the script.
   - `photometry/landolt.py:33` depends on
     `astrometry.services.propagate_proper_motion`.
   - The photometry QA uses `astrometry.qa.plot_solution`.
   - Mosaic-specific algorithms are in `pipelines/mosaic.py`: `_background`,
     `prepare_image`, `_source_mags`, `measure_offsets`.
4. **Reading a table can write to it.**
   - `load_setup` reads the table through `build_metadata`
     (`reduction.py:225`).
   - That path can rewrite the file's metadata: the `rdxdir` at
     `metadata.py:1148-1155`, and again at `reduction.py:239`.

   Reading should go through `read_metadata`, with any update explicit and
   logged.

### 3.3 Proposed layout

```
pipnick/
  pkg/        exceptions, logger, provenance.py (hashes), paths.py (product names)
  datamodel/  dataset.py (DataSet, Product), image.py (ReducedImage),
              catalog.py (SourceCatalog), psf.py (PSFModel),
              calibration.py (Calibration), settings.py (dataclasses)
  metadata/   table.py, setups.py, groups.py (comb_id, mosaics), build.py, names.py
  image/      io.py (one reader), background.py, detect.py
  reduction/  processing.py (frame processing, combination), combiner.py
  astrometry/ fit, match, services, cache.py, solve.py, header.py
  photometry/ moffat, psf_photometry, aperture, psfstars.py (was starfind),
              landolt, extinction, calibrate.py (fit), standards.py
  mosaic/     coadd, psfmatch, scaling, prepare.py
  utils/      crossmatch.py, lsq.py, masking, rgb
  qa/         display.py, sources.py, psf.py, calibration.py, mosaic.py
  pipelines/  base.py, reduction.py, astrometry.py, photometry.py,
              calibrate.py, mosaic.py   (drivers only, about 150-400 lines)
  scripts/    unchanged
```

Rules:

- Pipelines import `datamodel` and the algorithm packages, but **no pipeline
  imports another pipeline**.
- `datamodel` depends only on `pkg` and the algorithm modules (e.g.
  `wcs_from_header`), never on `pipelines`. This avoids circular imports.

### 3.4 Splitting the large modules

| Module | Split into |
|---|---|
| `pipelines/astrometry.py` | `astrometry/cache.py` (`_anet_*`, `_record_anet`, `_restore_solution`, `_invalidate_chains`, `:487-652`), `astrometry/header.py` (`_remove_wcs`, `_solution_header`, `apply_wcs`, `:249-369`), `astrometry/solve.py` (grouping and solving strategy, `:655-1341`), `qa/` (reporting); the driver stays (about 180 lines) |
| `pipelines/calibrate.py` | `datamodel/calibration.py` (TOML I/O, `:641-766`), `photometry/standards.py` (`match_standards`, `neighbors`, `_measurement_table`, `:202-397`), the band fit next to `photometry/calibrate.py`, the application (`apply_to_*`) in its own module, and QA in `qa/calibration.py` |
| `pipelines/reduction.py` | planning (`:243-659`) to `datamodel`; provenance (`:518-824`) to `pkg/provenance.py`; frame processing (`:893-1225`) to `reduction/processing.py` |
| `pipelines/mosaic.py` | preparation and scaling (`:336-516`) to `mosaic/prepare.py`; `check_mosaic_groups` to `metadata/groups.py` |
| `metadata.py` | as in the layout above; make `_safe_name`, `_filter_name` and `_window_name` public |

---

## 4. Redundancy and consolidation

### 4.1 Source detection: one function, used with different settings

Sources are detected in five places, with four implementations:

| # | Location | Method |
|---|---|---|
| 1 | `astrometry/detect.py:12` | global clipped statistics, then `DAOStarFinder`, brightest `nmax` (on data with masked pixels set to 0) |
| 2 | `photometry/starfind.py:101` | global statistics, then `DAOStarFinder` at `psf_thresh`, then Kron, isolation and flux cuts (the PSF stars) |
| 3 | `photometry/psf_photometry.py:205-211` | `Background2D` and `MADStdBackgroundRMS`, then `DAOStarFinder` inside `IterativePSFPhotometry` (also run on the residual images of later iterations) |
| 4 | `pipelines/mosaic.py:349-351` | `photutils.segmentation`, to mask sources for the constant background |
| 5 | `pipelines/mosaic.py:661` | `build_moffat_psf_model` again, on the matched image (implementation 2 reused) |

**How they differ.** The detection itself is the same in 1–3: `DAOStarFinder`
at a threshold times the noise, for a given FWHM. They differ in four ways:

- **Background and noise:** a global sigma-clipped median and standard
  deviation in 1 and 2, but a 2D background and a MAD noise estimate in 3.
- **Masked pixels:** set to 0 in 1; carried as a mask in 2 and 3.
- **Parameters:** threshold, minimum separation, and whether sources at the
  border are excluded.
- **Selection afterwards:** the brightest `nmax` in 1; Kron radii, isolation
  and flux cuts in 2; none in 3.

None of these differences is a reason for separate code; each is a parameter
or a selection.

**Proposal.** Write one function, and keep the selections as separate
functions that operate on its result.

- **`image/detect.py`:**
  `find_sources(data, mask=None, error=None, fwhm=..., thresh=5., background='2d', box_size=64, min_sep=1., exclude_border=False, max_sources=None) -> SourceCatalog`.
  - `background` is `'global'` (clipped median and standard deviation) or
    `'2d'` (`Background2D` and `MADStdBackgroundRMS`).
  - Masked pixels are always passed as a mask, never zeroed.
  - It returns the positions, fluxes, peaks and sharpness/roundness from the
    finder, plus the background and noise it used, in the catalogue's
    metadata. The catalogue is sorted by flux.
- **Each caller runs it as often as it needs, with its own settings:**

  | Caller | Settings and selection |
  |---|---|
  | Astrometry | `thresh` high, `exclude_border=True`, `max_sources` (replaces 1) |
  | PSF stars | `thresh=psf_thresh`, then `photometry/psfstars.py: select_psf_stars(catalog, data, fwhm, min_dist, min_flux)`, which keeps the Kron, isolation and flux cuts (replaces the detection in 2) |
  | Photometry | a low threshold for a complete census (replaces 3). `IterativePSFPhotometry` accepts as its *finder* any callable that takes `(image, mask=...)` and returns a table with x and y positions (e.g. `x_centroid`), so wrap `find_sources` in a small callable class (`image/detect.py: SourceFinder`) that holds the settings. The residual images of later iterations are then searched with the same function and settings. |
  | Mosaic QA | refit of the matched PSF (5): reuses the PSF-star path above |

- **Selections are functions or methods of `SourceCatalog`** (section 5.3):
  - `brightest(n)`;
  - `above_snr(x)`;
  - `isolated(min_dist)`;
  - `unflagged()`;
  - `inside(shape, margin)`.

  A step can cut one catalogue in different ways, or run the detection
  again with different settings, whichever is more natural; neither requires
  detecting only once.
- **The mosaic's background mask** (4) is a mask, not a catalogue, and
  `photutils.segmentation` is the right tool for it. Keep it, but move it to
  `image/background.py` (4.2) as the `mask_sources` option of the background
  estimate.
- **Optional reuse.** Because the photometry catalogue doesn't depend on the
  WCS (the celestial coordinates are added at the end,
  `pipelines/photometry.py:37`, and refreshed when the WCS changes,
  `:459-473`), a step could also take the brightest sources of an existing,
  current photometry table instead of detecting again. This is an option,
  not a requirement of the design.

**Risk: the astrometry.net cache.** Its key (`_input_hash`,
`astrometry.py:372`) includes the rounded positions of the detected sources.
If moving the astrometry to `find_sources` changed its detections (e.g., by
masking rather than zeroing the bad pixels), the keys would change, and an
unsolved image could be resubmitted.

- With `background='global'` and the astrometry's settings, `find_sources`
  should reproduce `astrometry/detect.py` exactly. Add a test that the two
  give identical catalogues before switching.
- If a behaviour change is wanted later, add a fallback that looks up the
  cache by image name.
- Under the standing rule, ask before any run that could submit a job.

### 4.2 Background estimation

There are four approaches:

- **Global clipped median:** `astrometry/detect.py:38`, `starfind.py:165`,
  `starfind.py:363`.
- **2D:** `psf_photometry.subtract_background` (`:101`, box 64), which the
  mosaic reuses with box 256 (`pipelines/mosaic.py:348`).
- **Constant with sources masked by segmentation:**
  `pipelines/mosaic.py:349-353`.
- **Local annulus, done four different ways:**
  - `kron_stats` (`:297`)
  - `extract_stamps` (`:388`)
  - `LocalBackground(MMMBackground)` (`psf_photometry.py:213`)
  - `ApertureStats.median` (`pipelines/mosaic.py:412-415`)

**Proposal:** `image/background.py` with two functions:

- `estimate_background(data, mask, method='2d'|'constant'|'none', box_size=64, mask_sources=True) -> (background, rms)`;
- `local_background(data, mask, x, y, r_in, r_out, statistic='median'|'mmm')`.

### 4.3 Cross-matching

| Function | Space | Matching |
|---|---|---|
| `astrometry/match.py:119 match_nearest` | pixels | KD-tree, nearest per reference |
| `photometry/landolt.py:283 match_landolt` | sky | nearest source per star, then unique per source |
| `mosaic/scaling.py:25 mutual_matches` | sky | mutual nearest neighbours |
| `pipelines/calibrate.py:1391 _cross_match` | sky | greedy linking into source IDs |

**Proposal:** `utils/crossmatch.py` with two functions:

- `match_sky(c1, c2, radius, unique='mutual'|'first'|'second') -> (i1, i2, sep)`;
- `link_sources(coords_list, radius) -> list of id arrays`.

`match_landolt` would become a proper-motion propagation followed by
`match_sky`. The pixel-space `match_nearest` stays separate.

### 4.4 Hashes and rerun logic

**Eleven hash helpers.**

- `reduction._table_hash`, `_raw_hash`
- `photometry._settings_hash`, `_image_hash`, `_wcs_hash`
- `astrometry._input_hash`
- `astrometry/fit.wcs_hash`
- `calibrate._hash`, `_file_hash`, `calibration_hash`
- `mosaic._hash`

Three of them are byte-identical: `photometry.py:86`, `calibrate.py:178`
and `mosaic.py:73`.

**Five mechanisms decide whether an output is up to date:**

| Step | Record | What is compared |
|---|---|---|
| Reduction | header keywords | `PNRAWHSH`, `PNTABH`, `PN<cal>H` |
| Photometry | table metadata | the image, settings and WCS hashes, with a separate coordinate-only update |
| Astrometry | an index file plus the header | the input hash, plus a chain check |
| Calibration | TOML | `input_hash`, plus `fit_hash` to detect edits |
| Mosaic | header keywords | `PNMOSHSH`, `PNGRPHSH` |

Also, `calibrate_all` changes `settings['input_hash']` on each night
(`calibrate.py:1215`).

**Proposal:**

- `pkg/provenance.py`, with `digest`, `json_hash`, `file_hash`, `names_hash`,
  `image_hash(header)` and `header_wcs_hash(header)`. The outputs must stay
  byte-identical, so existing products remain up to date; a test should
  enforce this.
- A `Provenance` record that each step builds (its input hashes and settings
  hash) and stores in a header or in table metadata.
- Each step answers `is_current(existing) -> (bool, reason)`. The astrometry
  chain and the edits to the calibration file are handled as overrides.

### 4.5 Finding and checking the input images

The same sequence appears four times:

1. `plan_products`, keep the science products;
2. if the image doesn't exist, warn and skip (`photometry.py:440`,
   `astrometry.py:450`, `calibrate.py:789`, `mosaic.py:154`);
3. if `PNSETUP` doesn't match the tag, raise (`photometry.py:445`,
   `astrometry.py:454`, `calibrate.py:796`, `mosaic.py:158`).

Other duplicated pieces:

- the check that the photometry table is current (`calibrate.py:799-809`,
  `mosaic.py:168-175`);
- the check that the reduction directory exists (`photometry.py:415`,
  `astrometry.py:1498`, `calibrate.py:1151`);
- two selection functions (`photometry._select:286`, which works on
  dictionaries, and `astrometry._select:945`, which sets attributes).

Consolidate these in `DataSet.science_images()`, `DataSet.select()` and
`ReducedImage.photometry()` (section 5).

### 4.6 Summaries, QA, input/output and small helpers

- **Summary writers.** Four `_summarize` functions:
  - `photometry.py:509`
  - `astrometry.py:866`
  - `calibrate.py:1563`
  - `mosaic.py:1007`

  Each builds a `Table`, converts `-1` sentinels to masked values, writes
  `summary.ecsv` and logs it. Replace them with a `Summary` class in
  `pipelines/base.py`, where `None` means masked. That also retires the `-1`
  sentinels still written to the mosaic `INPUTS` table
  (`pipelines/mosaic.py:624-628`).
- **QA.**
  - Three images are displayed with ZScale: `astrometry/qa.py:43`,
    `pipelines/mosaic.py:1034` and `photometry/fit.py:153`.
  - The same save-a-figure sequence (`Figure`, `tight_layout`, `mkdir`,
    `savefig`, log) appears in `calibrate.py:1615-1660` and
    `mosaic.py:1032-1061`.

  Consolidate them in `qa/display.py`: `show_image(ax, data, mask)` and
  `save_figure(fig, ofile)`.
- **Reading reduced images.** Six readers, which differ in how they treat
  the mask and the error:
  - `photometry.measure_image` (`CCDData.read` plus `fits.getheader`, so the
    file is opened twice)
  - `astrometry._read_image`
  - `mosaic.prepare_image`
  - `utils/rgb.read_image`
  - `starfind.as_masked_array`
  - `reduction._calibration_data`

  The test that the error is finite and positive is written twice
  (`photometry.py:147`, `mosaic.py:376`). Proposal:
  `image/io.read_reduced(path) -> (data, mask, error, header)`.
- **Least-squares fits.** `fit_calibration` (`photometry/calibrate.py:234`),
  `residual_drift` (`:358`) and `solve_offsets` (`mosaic/scaling.py:126`)
  each build and solve the normal equations, and each scales the covariance
  by the reduced χ². Proposal:
  `utils/lsq.weighted_lstsq(X, y, err, gpm, scale_cov=True)`. The rejection
  policies differ, so they stay with the callers.
- **The Moffat FWHM** is computed in several places:
  - `FitMoffat2D.to_fwhm`, duplicated in `FitEllipticalMoffat2D`;
  - `psf_photometry.gamma_to_fwhm` and `process_par`;
  - `EllipticalMoffat2DModel.fwhm` and `EllipticalMoffatPSF.fwhm`;
  - `mosaic/psfmatch.moffat_fwhm`.

  **This is an inconsistency, not just duplication.** For an elliptical PSF,
  `process_par` uses the *mean* of the axes (`psf_photometry.py:277`), but
  `moffat_fwhm` uses the *major* axis (`psfmatch.py:66`). So `fwhm_in` in
  the mosaic disagrees with `fwhm_pix` in the photometry table. Resolve it in
  a `PSFModel` class (5.4) with one documented definition, and give the major
  and minor axes explicitly.
- **Names of products** are built in several places:
  - `{stem}_phot.ecsv` in `photometry.py:434`, `calibrate.py:788` and
    `mosaic.py:164`;
  - the astrometry sidecar files, built differently at `astrometry.py:365`
    and `:593`;
  - the calibration files.

  Proposal: `pkg/paths.py`.
- **Small duplicates:**
  - `_night` (`calibrate.py:188`), inlined at `mosaic.py:194`;
  - an MJD that is masked becoming `None` (`photometry.py:477`,
    `mosaic.py:181`);
  - the pixel scale (`photometry.py:449`, `astrometry.py:463`);
  - the pointing as a `SkyCoord` (`astrometry.py:460`, `reduction.py:877`);
  - the overlap of two bounding boxes (`mosaic._check_overlap`,
    `coadd.overlap_differences`).

---

## 5. Object orientation

Data currently moves between functions as:

- a five-element tuple from `load_setup`, unpacked identically in five
  places;
- product dictionaries from `plan_products`;
- image dictionaries in the calibration and mosaic steps;
- settings dictionaries.

Only the astrometry uses classes (`AstrometryImage`, `_Session`).

The proposed classes wrap these structures. They keep the on-disk formats
unchanged, and they wrap `Table` objects rather than subclassing them.

### 5.1 `Product` and `DataSet`

**`Product`** (a frozen dataclass):

- Replaces the dictionary built in `plan_products`.
- Fields: `frametype`, `comb_id`, `indx`, `filter`, `exptime`, `file`.

**`DataSet`:**

- Replaces the `load_setup` tuple.
- Attributes: `metafile`, `table`, `rdxdir`, `tag`, `camera`.
- Methods:

  | Method | Absorbs |
  |---|---|
  | `from_file(table, rdxdir=None)` | `load_setup`, `check_single_setup`, `check_comb_groups`, `check_mosaic_groups` |
  | `products(frametype=None)` | `plan_products`, `_group_indices` |
  | `table_values(product)` | `table_values` |
  | `raw_files(product)` | — |
  | `pixel_scale(product)` | — |
  | `pointing(product)` | — |
  | `night(product)` | — |
  | `select(items)` | both `_select` functions |
  | `outdir(step)` | — |
  | `science_images()` | (section 4.5) |

- Keep `load_setup` and `plan_products` as thin wrappers until the tests
  move over.

### 5.2 `ReducedImage`

- Attributes: `path`, `product`, `dataset`, and the header, loaded when
  first needed.
- Methods:
  - `read()` returns the data, mask and error (section 4.6).
  - `check_setup()`.
  - `photometry(require=True)` returns the `SourceCatalog`, after checking
    that it is current. It absorbs `calibrate._collect_images` and
    `mosaic.collect_images`.
- Properties: `wcs`, `has_solution`, `image_hash`, `wcs_hash`, `phot_file`,
  `zero_point`.

`AstrometryImage` would hold a `ReducedImage` plus its solving state, and the
image dictionaries in the calibration and mosaic steps would become
`ReducedImage` objects.

### 5.3 `SourceCatalog`

This wraps a `_phot.ecsv` table.

- Properties: `x`, `y`, `coords`, `psf`.
- Methods:

  | Method | Absorbs |
  |---|---|
  | `from_measurement` | the table assembly in `measure_image` |
  | `update_coordinates(header)` | `celestial_coordinates` |
  | `instrumental_mag(kind)` | `calibrate._instrumental`, `photometry.py:172-180`, `mosaic._source_mags` |
  | `usable()` | — |
  | `neighbors(i)` | `calibrate.neighbors` |
  | `brightest(n)`, `above_snr(x)`, `isolated(d)`, `unflagged()`, `inside(shape, margin)` | the selections now in `astrometry/detect.py` and `starfind.py` (section 4.1) |
  | `is_current(image)` | — |
  | `read`, `write` | — |

### 5.4 `PSFModel`

- Attributes: `fittype`, `par`.
- Constructors: `from_fit`, `from_meta`, `circular(fwhm, beta)`.
- Properties: `fwhm` (one documented definition), `fwhm_major`,
  `fwhm_minor`, `beta`.
- Methods: `to_photutils()`, `render(size, scale)`, `to_meta()`.
- It absorbs:
  - `gamma_to_fwhm`
  - `process_par`
  - `moffat_fwhm`, `moffat_beta`, `render_psf`, `target_psf`
  - `mosaic._psf`
  - the photutils model construction in `psf_analysis:216-221`
  - the duplicated `to_fwhm` and `to_gamma` static methods

Choosing one FWHM definition changes the mosaic target, so the affected
mosaics will be rebuilt; note this in the release notes.

### 5.5 `Calibration`

This wraps the calibration TOML.

- Attributes: `path`, `content`, `edited`.
- Methods:
  - `read` and `write` absorb `read_calibration`, `write_calibration`,
    `_toml_text`, `_toml_value` and `calibration_hash`.
  - `bands()` returns `BandCalibration` objects (reusing `CalibrationFit`),
    and absorbs `_calibrations` and `_extinction_source`.
  - `zero_point(band, airmass)`, `summary_rows()`, `is_current(input_hash)`.

### 5.6 Settings dataclasses

`PhotometrySettings`, `AstrometrySettings` (replacing the parameters held by
`_Session`), `CalibrationSettings` and `MosaicSettings`. Each is frozen, with:

- `validate()`, which absorbs the option checks in `mosaic.py:881-889` and
  `calibrate.py:1153-1156`;
- `hash()`, which is the `json_hash` of the dataclass.

Making them frozen also removes the in-place change of the settings at
`calibrate.py:1215`.

### 5.7 `PipelineStep`

A base class in `pipelines/base.py`:

- Class attributes: `name` (also the output subdirectory) and
  `summary_columns`.
- `__init__(dataset, settings, overwrite, qa)`.
- Methods: `inputs()`, `is_current(item)`, `process(item) -> row`,
  `write_qa(item, result)`.
- A template method, `run(selection=None)`, provides the common loop:
  1. find the inputs;
  2. skip any that are current;
  3. process the rest, catching and recording errors (which also fixes
     section 6.2, item 6);
  4. write the summary, the verdict and the next-step hint;
  5. write the QA figures.

The steps fit it as follows:

- **Reduction and photometry:** one item per image.
- **Calibration:** one item per night.
- **Mosaic:** one item per mosaic, with the WCS of each group computed
  first. Its `inputs()` takes several `DataSet`s.
- **Astrometry:** its solving order is global across pointing groups, so it
  overrides `run`, but it reuses `inputs`, `Summary` and the QA.

### 5.8 Migration order and risks

1. `pkg/provenance.py`, `pkg/paths.py`, `utils/lsq.py`,
   `utils/crossmatch.py`. These are pure refactoring, guarded by tests that
   the hashes are byte-identical.
2. `Product` and `DataSet`, keeping wrappers.
3. `read_reduced` and `ReducedImage`, then the functions that find the
   inputs.
4. `PSFModel`, once the FWHM definition is decided.
5. `SourceCatalog`, `Summary`, and `qa/`.
6. Split the large modules (3.4).
7. The common detection function (4.1). Moving the astrometry to it must not
   change its detections, or the astrometry.net cache is invalidated.
8. `PipelineStep`.

**Risks:**

- The header keywords and the ECSV and TOML metadata are the on-disk
  contract, and must not change.
- The tests use private names. Update them as the names become public.
- Circular imports. `datamodel` must not import `pipelines`.

---

## 6. Log messages

### 6.1 Current state, by step

**Logging infrastructure** (`scripts/scriptbase.py`, `pkg/logger.py`):

- The console default is INFO.
- **No log file is written unless `--log_file` is given**
  (`scriptbase.py:103-106, 130-135`). The default name goes in the current
  directory, not the reduction directory.
- **No step records the pipnick version.**
- `pipnick_metadata` reports argument errors and then returns, so the exit
  status is 0 (`scripts/metadata.py:92-115`). Scripts can't detect the
  failure.

**`pipnick_metadata`** is good:

- `summarize` (`metadata.py:925-982`) prints the setups, the counts by frame
  type, every commented-out row with its reason, and the missing
  calibrations.

Gaps:

- It gives no next-step hint.
- It doesn't say that the checks that need the pixel data (bias level,
  saturation) were skipped without `--classify stats`.
- The warning about a changed `rdxdir` is split across lines
  (`:1152-1155`).

**`pipnick_reduce`** is the weakest step:

- **The run header** (`reduction.py:126-129`) gives the file, the setup and
  the output directory, but no frame counts by type, no list of the products
  to build, and nothing on which are up to date.
- **Each raw frame produces 5–7 INFO lines** (`:922`, `:1037`, `:1040`,
  `:1048`, `:1056`, `:1064`, and the clipping lines). Most don't name the
  frame, and about 400 lines a night bury the warnings.
- **There is no progress count.** The cosmic-ray step is slow.
- **Measurements are reported without a verdict.** The bias level and read
  noise (`:1202, 1221`) aren't compared with the camera's `[quality]` range
  or its nominal read noise. Nothing at all is reported for the flats (level,
  illumination, saturation) or for the science frames (sky, saturated
  fraction, cosmic-ray fraction).
- **There is no closing summary.** Problems that degrade the output are only
  warned about mid-stream: no bias (`:151`), no flat (`:848`), unknown read
  noise (`:1042`), offsets combined without registration (`:886`).
- **It writes no summary table and no QA figure.**
- **One failing frame stops the whole run.**
- **What works well:** the reasons for reprocessing (`:736-756, 1274`) and
  the provenance written to the headers.

**`pipnick_astrometry`** is the clearest step:

- each failure has a code;
- `summary.ecsv` gives a reason and a suggestion for every unsolved image;
- the warning about a pointing error says what to fix and which script to
  rerun.

Gaps:

- The RMS has no context: is 0.42″ good?
- There is no progress count.
- "Rejected the solution …" is logged at INFO (`_fail`, `:181`), even when it
  decides the outcome.
- QA is written only for images newly solved in this run (`:927-929`), so
  unsolved images, which most need a look, get none.

**`pipnick_photometry`:**

- The per-image line (`photometry.py:219-225`) and the summary are useful.
- The lower-level messages are jargon, at INFO, without the image name:
  - "Median and standard deviation: 412.3 9.87" (`starfind.py:166`, no
    units);
  - "Mean Kron radius of sources: 3.4";
  - "Consolidating groups [3, 7, …]" (`psf_photometry.py:307`), which can be
    long.
- The FWHM is logged twice.
- The sources removed by `vet_sources` are only counted (`:279`), then
  dropped, so the user can't see which detections were removed.
- Only `ValueError` is caught (`:482`).
- QA is written only for images that are measured again (`:492`).

**`pipnick_calibrate`** is well instrumented:

- every standard measurement has a reason;
- failures and the extinction fallback are warned about, with a banner at
  the end;
- the drift is tested for significance;
- an edited calibration file is reported when it is kept.

Gaps:

- The summary has no zero-point column, so a night 0.5 mag low (clouds)
  isn't visible there.
- The summary's `file` column holds the action taken (fit, up to date,
  edited and kept).
- The fitted extinction isn't compared with the nominal value.

**`pipnick_mosaic`:**

- The per-mosaic line, the seeing-outlier warning and the `PNWARN` header
  cards are useful.
- **The `failed` status can't actually occur.** Errors aren't caught, and
  `build_mosaic` raises when an image doesn't overlap (`:582`), so one bad
  image stops every mosaic and no summary is written.
- The warning for the Gaussian fallback doesn't name the image
  (`psfmatch.py:290-293`).
- The scale and FWHM of each image are in the `INPUTS` table, but never
  logged.

**Across steps:**

- The status words differ between steps:

  | Step | Status values |
  |---|---|
  | Astrometry | `solved`, `unsolved`, `pending`, `skipped` |
  | Photometry | `measured`, `up to date`, `coordinates updated`, `missing`, `failed` |
  | Calibration | `ok`, `failed`, plus the `file` column |
  | Mosaic | `constructed`, `up to date`, `failed` (never set) |
  | Reduction | none |

- The paths are written in different ways: "Saved", "Wrote", absolute,
  relative and bare names.
- No step gives a next-step hint.
- QA is off by default everywhere.

### 6.2 Proposed messages

The aim is for an inexperienced observer to be able to answer four
questions from the log alone: what was done, was it OK, what needs
attention, and what do I run next.

| Step | Location | Level | Proposed message | Why |
|---|---|---|---|---|
| all | `ScriptBase.init_log` (new helper) | INFO | `pipnick 0.1.0 \| pipnick_reduce \| 2026-10-26 04:12 UT \| log: <rdxdir>/logs/pipnick_reduce_20261026-0412.log` | Version and log location; write a log by default. |
| all | each `*_all`, at the start | INFO | `Inputs: <table> (setup 2x2_slow_full, NickelC2); output: <dir>` | A consistent run header. |
| metadata | end of `build_metadata` | INFO | `Next: review n1_2x2_slow_full.ecsv (comment out bad frames with '#', check the object names and comb_id groups), then run: pipnick_reduce -t n1_2x2_slow_full.ecsv` | Next step. |
| metadata | `build_camera_datatable` (`classify='header'`) | INFO | `Checks that need the pixel data (bias level, flat illumination, saturation) were skipped; use --classify stats to apply them.` | Says what wasn't checked. |
| reduce | after `plan_products` | INFO | `Setup A: 9 bias, 10 sky flats (B, V, R), 22 science frames in 14 groups → 17 products: 3 up to date, 14 to build (12 new, 2 because the bias changed).` | The plan up front. |
| reduce | each product | INFO | `[5/17] Building PG1530+057_V_C014_2024-06-25.fits (science, V, 3 frames, 60 s) ... done in 41 s` | Progress and timing. |
| reduce | `process_frame` | DEBUG (from INFO) | `d1047.fits: overscan 301.4 ADU; trimmed to 1024x1024; bias and SkyFlat_V applied; 812 cosmic-ray pixels (0.08%)` | One line per frame, naming the frame; the intermediate lines move to DEBUG. |
| reduce | `_report_bias_stats` | INFO / WARNING | `Bias: level 301–304 ADU (expected 280–400) OK; read noise 2.74 e- (nominal 2.70) OK.` or `Read noise 7.9 e- is 2.9× the nominal 2.70 e- for the slow readout; check for pattern noise or a mislabeled readout speed before trusting the errors.` | A verdict against the camera's values. |
| reduce | flats | INFO / WARNING | `SkyFlat_V: 10 frames, levels 18,200–31,500 e- (scaled to 24,100); normalized flat 0.91–1.06 in the central 90%.` / `d1012 has 9% saturated pixels` | The flats aren't checked at all now. |
| reduce | science | INFO | `PG1530+057_V_C014: sky 412 e-/pix (6.9 e-/s), 0.02% saturated, 812 cosmic-ray pixels.` | The sky level shows twilight, moonlight or clouds. |
| reduce | `_select_flat` | WARNING | `No V flat in this setup: the V science images (C014, C015) are NOT flat-fielded, and their photometry will have spatial errors of ~5–10%.  Add V flats (or flats from another night) to the table and rerun.` | Consequence and remedy. |
| reduce | end | INFO / WARNING | `REDUCTION OK: 17 products (3 up to date, 14 built); summary in reduce/summary.ecsv.  Next: pipnick_astrometry -t <table>` or `REDUCTION FINISHED WITH 2 PROBLEMS: no B flat (4 images not flat-fielded); 1 group combined with offsets.  See the warnings above.` | Closing verdict, repeating the problems. |
| reduce | each product | WARNING | `Failed to build <file>: <error>.  Continuing with the other products; see reduce/summary.ecsv.` | One failure no longer stops the run. |
| astrometry | solved | INFO | `[7/22] Solved <image> from Gaia: 48 stars, RMS 0.21 arcsec (typical 0.1–0.4) OK` | Progress and context. |
| astrometry | `_fail` (final rejection) | WARNING | (current text, raised to WARNING when it decides the outcome) | The decisive rejection becomes visible. |
| astrometry | end | INFO / WARNING | `ASTROMETRY OK: 22/22 solved.  Next: pipnick_photometry -t <table>` / `NEEDS ATTENTION: 2 unsolved (see the suggestions above; QA in astrometry/qa)` | Verdict and next step. |
| photometry | `starfind.detect_sources` | DEBUG | `<image>: background 412.3 e-, noise 9.87 e- (clipped)` | Units and the image name; moved to DEBUG. |
| photometry | `consolidate_groups` | DEBUG | `Merged 14 groups of blended sources` | No long id lists at INFO. |
| photometry | `measure_image` | WARNING | `<image>: PSF FWHM 4.1 arcsec is well above typical Nickel seeing (1.5–2.5 arcsec); check the focus or tracking.` / `only 3 stars used to model the PSF; the PSF magnitudes may be unreliable.` | Sanity checks; the expected values come from the camera. |
| photometry | `vet_sources` | INFO | `Removed 7 of 312 sources (5 moved >1 FWHM, 2 non-positive flux); listed in photometry/<image>_removed.ecsv` | Lets the user trace removed sources. |
| photometry | end | INFO | `PHOTOMETRY OK: 22 images (median FWHM 1.7 arcsec, 40–250 sources).  Next: pipnick_calibrate -t <table> (standards) or pipnick_mosaic` | Verdict and next step. |
| calibrate | per band | INFO / WARNING | `V: zp 24.31 ± 0.01 (typical 24.2–24.5), k = 0.14 mag/airmass (nominal 0.18), 34 stars, rms 0.023 mag OK` / `V zero point 0.6 mag below typical: likely clouds; consider not calibrating this night.` | Sanity check against expected values. |
| mosaic | each image | INFO | `[2/6] <image>: scale ×1.032 ± 0.004 (41 stars), FWHM 1.8 → 2.4 arcsec (Moffat kernel)` | Per-image result. |
| mosaic | `image_kernel` fallback | WARNING | Name the image in the current message. | Can't be attributed now. |
| mosaic | each mosaic | WARNING | `Mosaic NGC_7331_V failed: <image> does not overlap the others; continuing with the other mosaics.` | Makes the `failed` status real. |

Most of the checks against expected ranges ("typical", "OK") need those
ranges in one place. Proposal: an `[expected]` section in the camera TOML,
with the read noise, sky rates by filter, FWHM range, zero point by band and
astrometric RMS (7.2, P6).

Also standardize the following:

- one verb for written files, e.g. `Wrote <path relative to rdxdir>`;
- the frame or image name in every per-image message;
- a final line in every step, beginning `<STEP> OK`, `<STEP> NEEDS ATTENTION`
  or `<STEP> FAILED`.

---

## 7. Evaluating and interrogating the results

### 7.1 Current state

What exists:

- **Summary tables:** astrometry, photometry, calibration and mosaic. None
  for the reduction or the metadata.
- **QA figures, all opt-in:**
  - astrometry: the solution vectors;
  - photometry: the PSF and the sources;
  - calibration: the residuals;
  - mosaic: a four-panel figure.
- **Header provenance:** the `PN*` keywords and `HISTORY` cards.
- **Calibration files:** the editable TOML, and the residual table, which
  gives a reason for each standard measurement.
- **Astrometry records:** the astrometry.net submission records.
- **The mosaic `INPUTS` table.**

Gaps:

1. There is no single view of a night across the steps.
2. QA is off by default. Even with `--qa`, nothing is drawn for products that
   are up to date or for images that aren't solved.
3. The reduction step writes no summary and no QA figure.
4. The summaries can't be joined:
   - the astrometry and photometry are keyed by image, the calibration by
     night and band, and the mosaic by mosaic;
   - the status words differ.
5. No pipnick version is recorded in the outputs.
6. Logs are lost by default.
7. A source can't be traced through the steps:
   - the photometry `id`s are per image;
   - the cross-image `source_id` exists only in the `_cal.ecsv` files;
   - the removed sources are discarded.
8. `docs/basics.rst` and `docs/example.rst` are out of date.

### 7.2 Proposals, by priority

**P1. `pipnick_status -t <table>` (read-only, fast).**

- One row per science image, joining the summaries of each step and the
  header keywords.
- Columns: `image, object, filter, exptime, airmass, reduced, wcs (method/rms), nsrc, fwhm, zp, calibrated, mosaic, attention`.
- Each step's cell uses one vocabulary: `ok | stale | missing | failed | n/a`.
  `stale` comes from the existing hash checks.
- It ends with the next command to run.

It uses only information that already exists, and answers "where am I, and
what's broken?"

**P2. Consistent outputs from every step.**

- QA on by default (with `--no_qa`), written for any product without a
  figure, including unsolved images.
- Each step writes `<step>/summary.ecsv`, with the common columns
  `image | status | reason | suggestion`, plus the closing verdict and the
  next-step hint (6.2).
- A new `reduce/summary.ecsv` with the bias level, read noise, flat levels,
  sky level, cosmic-ray fraction and the calibrations used.
- A reduction QA page: the bias, the flat for each filter, and a thumbnail of
  each science image.

**P3. Records of each run.**

- Logs written by default to `<rdxdir>/logs/`.
- `PNVERS` (the pipnick version) in every product header and in the
  metadata of every summary.
- A `<rdxdir>/logs/runs.jsonl` file, with one record per run: the script,
  its arguments, the version, the start and end times, the counts by status,
  and the warnings.
- A non-zero exit status on errors.

**P4. `pipnick_report`.**

- Renders P1 and the summaries into a static `<rdxdir>/report/index.html`.
- Contents:
  - an overview table of the night, colour-coded by status;
  - a page for each image with its QA figures and provenance, from the raw
    frames through the calibrations, WCS, PSF and zero point to the mosaic;
  - the calibration residuals, and the extinction banner when it applies;
  - links to the summary tables and the log.

This is the most useful single addition for a new observer, and for anyone
reviewing someone else's run.

**P5. `pipnick_inspect <image>`.**

- Prints, for one product:
  - its provenance: setup, `comb_id`, raw frames, the bias and flat, and
    whether each is up to date;
  - the WCS: method, reference, matches, RMS, pointing error;
  - the PSF: FWHM, number of stars, parameters;
  - the photometry counts, including the removed sources;
  - the calibration keywords and file;
  - the mosaics it belongs to, with its scale and FWHM in each.
- `--source RA,Dec` (or `--id`) traces one star through the `_phot.ecsv`,
  `_cal.ecsv` and `_residuals.ecsv` files. For a standard star, it shows
  whether it was used, rejected or excluded, and why.
- `--display` shows the image with its sources overlaid. The existing
  display code (1.1) could be revived here, written against the new data
  model.

**P6. Expected ranges in the camera TOML (`[expected]`).** These feed the
"OK / needs attention" checks in 6.2, so the definition of normal for each
camera lives in one place.

**P7. Small fixes.**

- Keep the removed sources in the photometry table with a `removed` flag,
  or in a separate file, rather than dropping them.
- Add a zero-point column to the calibration summary, and rename its `file`
  column to `action`.
- Mask the missing values in the mosaic `INPUTS` table instead of writing
  `-1`.
- Rewrite `docs/basics.rst` as a "reading the logs and summaries" page, with
  a table of typical values for the Nickel cameras.

---

## 8. Suggested order of work

| Step | Work | Effect | Risk | Status |
|---|---|---|---|---|
| 1 | Remove the unused code and dependencies (§1); fix the convention items (§2), including making the shared private helpers public | cleaner code; no `loess` packaging problem | low | **complete** |
| 2 | `pkg/provenance.py`, `pkg/paths.py`, `utils/lsq.py`, `utils/crossmatch.py` (§4) | less duplication; one hash implementation | low (needs byte-identical hash tests) | **complete** |
| 3 | `DataSet`, `Product`, `ReducedImage`, `read_reduced`; `PSFModel` with one FWHM definition (§5) | removes the coupling between pipelines and the FWHM inconsistency | medium (some mosaics are rebuilt) | **complete** |
| 4 | Logging (run header, progress, verdicts, next steps, default log files, version); a reduction summary; `pipnick_status` (§6, §7 P1–P3) | the most direct benefit to new observers | low | **complete** |
| 5 | Split the large modules; `SourceCatalog`, `Summary`, `qa/`; QA on by default; `pipnick_report`, `pipnick_inspect` (§3, §5, §7) | maintainability and interrogation | medium | **complete** |
| 6 | The common detection function and catalogue selections (§4.1); `PipelineStep` (§5.7) | one detection implementation, used with different settings by every step; a common step loop | higher (the astrometry.net cache, unless the astrometry's detections are reproduced exactly) | |

### Progress log

#### Step 1: unused code and conventions (complete, 2026-10-26, branch `review`)

Checks after the work:

- `tox -e 3.13-test,codestyle` passes in a fresh environment: 161 tests
  passed, 1 skipped.
- The docs build without warnings.
- All eight scripts start.

**Unused code removed (§1).**

- Modules: `utils/graphs.py`, `utils/log.py`, `utils/display_fits.py` and
  `utils/fits_class.py`, with their API pages.
- `photometry/fit.py`: `get_param_list`, `get_apertures` and `plot_sources`;
  `psf_plot` is kept.
- Members and helpers: `FitEllipticalMoffat2D.get_orig_phi`,
  `pkg.logger.clear_text_color` and `PipnickCodingError`.
- Dead code:
  - the commented-out code in `utils/masking.py`, `scripts/__init__.py` and
    `scripts/reduce.py` (the `--display` option and its import);
  - two commented-out lines in `pipelines/reduction.py`;
  - the unused `IPython` import in `tests/test_mask.py`.
- Unused imports and variables left behind (found with `pyflakes`), and an
  f-string without placeholders.
- `photometry/landolt.py: build_catalog` is kept, as a maintenance tool.

**Dependencies.** `IPython`, `loess`, `pandas` and `requests` were removed
from `pyproject.toml`. `licenses/README.rst` is updated: all remaining
dependencies are permissive, and a warning against restrictive dependencies
is kept.

**Conventions fixed (§2).**

- **Private names across modules.** The shared helpers are now public:
  - `metadata.safe_name` and `filter_name`;
  - `pipelines.photometry.image_hash` and `select_products`;
  - `solution_hash`, which was `_wcs_hash`. It is renamed to avoid confusion
    with `astrometry.fit.wcs_hash`, which hashes a WCS object.

  Each has a numpydoc docstring, and the metadata functions are in
  `metadata.__all__`.
- **Wrapping style.**
  - All hanging indents and trailing commas, in the package and the tests.
  - The one-element-per-line collections and signatures (e.g.
    `metadata_columns`, `logger.py`, `exceptions.__all__`, and the code adapted
    from astropy in `moffat.py`).
  - The backslash continuation in `services.py`.
- **Single fixes:**
  - the truthiness check (`pipelines/mosaic.py`);
  - the line over 99 characters (`moffat.py`), now a shorter URL;
  - the late imports in `pipnick/__init__.py`;
  - the `import ... as` aliases in `test_calibrate.py` and `test_mosaic.py`;
  - the docstring roles and package-name literals;
  - the unsorted glob in `test_astrometry.py`;
  - a module docstring for `utils/masking.py`, with the attribution of its
    polygon functions to `mangadap` (`mangadap/util/geometry.py`, BSD
    3-Clause).  The author confirmed `mangadap` as the origin.  The license
    is in `licenses/MANGADAP_LICENSE.rst`, listed in `licenses/README.rst` and
    in the `license-files` of `pyproject.toml`.
- **Updated:** the release notes.

**Changes beyond the review.**

- `pipnick_version --dependencies` also reports `astroquery`, `reproject` and
  `speclite`.
- **The `examples/` directory was removed by the author.** It held `all.py`,
  the notebooks `all.ipynb`, `reduction.ipynb`, `astrometry.ipynb`,
  `photometry.ipynb` and `final_calib.ipynb`, and 19 raw C2 frames in
  `data_example/raw/` (about 41 MB). Before the removal, all of them were
  confirmed to be defunct:
  - they imported modules that no longer exist (`pipelines.final_calib`,
    `utils.log`, `utils.display_fits`, `pipnick.convenience`), and
    `photometry.fit.plot_sources`, which was removed;
  - they called `reduce_all`, `astrometry_all` and `photometry_all` with
    arguments that were removed (`excl_files`, `excl_filts`, `excl_objs`,
    `use_table`, `save_inters`, `plot_final`, `plot_inters`), or with a
    directory where a metadata table is now required.

  Nothing else used the example data: the tests are synthetic.

  The `prune examples` line was removed from `MANIFEST.in`. Its
  `global-exclude` of `*.ipynb` files is **kept on purpose**: notebooks remain
  excluded from the source distribution for the time being.

**Left as is:**

- **Four `if` conditions** in `pipelines/mosaic.py`, `pipelines/astrometry.py`
  and `cameras/cameras.py` open a call or a parenthesis at the end of the
  line. Their contents are on the lines below, which arguably meets the
  convention (§2.2, uncertain).
- **Other one-element-per-line constructs** flagged by the automatic check
  (e.g. long `help` strings) have elements too long to pack.
- **The stale usage documents** `docs/basics.rst` and `docs/example.rst`
  still refer to removed code. They already referred to modules that no
  longer exist (`pipnick.convenience`), and their rewrite is deferred.
  `docs/example.rst` also refers to the example data removed with
  `examples/`.
- **The tests** still use private names in 30 places. The new rule exempts
  tests.

#### Step 2: pure helpers (complete, 2026-10-27, branch `consolidate_helpers`)

Checks after the work:

- `tox -r -e 3.13-test,codestyle` passes in a fresh environment: 176 tests
  passed, 1 skipped.
- The docs build without warnings.
- All eight scripts start.
- **The hashes are unchanged.**
  - The new tests compare every hash with reference values computed by the
    implementations that were replaced.
  - On a copy of the 2026-08-17 n3 (fast) reduction, `reduce_all` and
    `photometry_all` report all 24 reduction products and all 19 source
    tables as up to date.

**`pkg/provenance.py` (§4.4).**

- It provides `digest`, `json_hash`, `file_hash`, `names_hash`,
  `values_hash`, `wcs_hash`, `image_hash` and `solution_hash`.
- These replace:
  - `reduction._table_hash` and `_raw_hash`;
  - `photometry._settings_hash`, `image_hash` and `solution_hash`;
  - `calibrate._hash` and `_file_hash`;
  - `mosaic._hash`;
  - `astrometry.fit.wcs_hash`.
- `astrometry._input_hash` (the astrometry.net cache key) and
  `calibrate.calibration_hash` keep their roles, built on `digest` and
  `json_hash`.
- **Deferred:** the `Provenance` record and `is_current` in the second half
  of §4.4.

**`pkg/paths.py` (§4.6).**

- It provides `step_dir`, plus `photometry_table`, `points_file`,
  `calibrated_table`, `calibration_file` and `residuals_file`.
- The products derived from one image are named from the image's path, so
  `calibrate._collect_images` and `apply_to_catalogs` no longer take an
  output directory.
- **Left as is:** the names used only within one step (the astrometry.net
  cache files, the QA figures, the summaries and the mosaic products).

**`utils/lsq.py` (§4.6).**

- `weighted_lstsq` solves the normal equations for `fit_calibration`,
  `residual_drift` and `solve_offsets`. Each caller keeps its own rejection
  and its own rule for scaling the covariance by the reduced χ².
  - `solve_offsets` scales after its rejection loop, using the final
    selection. When the loop reaches `maxiter`, that selection includes one
    rejection that the solution doesn't. This behaviour is kept unchanged.
- On 300 random problems, the results agree with the old code to better
  than 1e-11 (relative). The differences are rounding from the reordered
  matrix products. The integer and boolean outputs are identical.

**`utils/crossmatch.py` (§4.3).**

- `match_sky(coords1, coords2, radius, unique='mutual'|'second'|'none')`
  serves three callers:
  - `match_landolt` (`'second'`);
  - the mosaic flux scaling (`'mutual'`; `mosaic.scaling.mutual_matches` is
    removed);
  - the linking below (`'none'`).
- `link_sources(ra, dec, radius)` replaces `calibrate._cross_match`. A list
  without celestial coordinates is passed fully masked.
- On 400 random crowded fields, all three give results identical to the old
  code.
- The pixel-space `astrometry.match.match_nearest` is unchanged, as
  proposed.

**Tests.** There are four new modules: `test_provenance.py`, `test_paths.py`,
`test_lsq.py` and `test_crossmatch.py`. `test_mutual_matches` moved to
`test_crossmatch.py`.

**Updated:** the release notes and the API pages.

#### Step 3: common objects (complete, 2026-10-28, branch `common_objects`)

Decisions made with the author before starting:

- **FWHM:** `PSFModel.fwhm` is the mean of the two axes, as the source
  tables already record it. The PSF matching of the mosaics explicitly uses
  `fwhm_major`, since it must broaden the widest axis, and `fwhm_in` is
  documented as the major-axis FWHM. Because both behaviours were kept, no
  product changed and no mosaic was rebuilt, contrary to the review's
  estimate (§5.4).
- **No wrappers:** `load_setup`, `plan_products`, `table_values` and
  `select_products` are removed, and every caller and test uses `DataSet`.
- **Scope:** `ReducedImage` is used by all four later steps.
- **Unchanged:** a reduction directory given to `DataSet.from_file` is still
  saved to the table (§3.2, item 4); to be revisited in Step 4.

Checks after the work:

- `tox -r -e 3.13-test,codestyle` passes in a fresh environment: 180 tests
  passed, 1 skipped.
- The docs build without warnings.
- All eight scripts start.
- **Readers:** on eight real reduced images (n2, n3), `read_reduced` returns
  the same data, mask, error and header as each of the four readers it
  replaced. The astrometry's inputs, and so its astrometry.net cache keys,
  are unchanged.
- **End to end:** on scratch copies of the 2026-08-17 n3 (fast) reduction,
  with the M11 images labelled as four mosaics in one WCS group, `main` and
  this branch ran the reduction, photometry, calibration (zero point,
  applied to the headers and catalogues) and mosaics.
  - All image, mask, error and table data are identical.
  - The only header and metadata differences are the calibration file's
    creation time (and so its hash) and hashes that include the absolute
    paths of the two copies.
  - Run on the copy built by `main`, this branch reports every reduction
    product, source table, calibration, header calibration and mosaic as up
    to date, and rebuilds nothing.

**`pipnick/datamodel/` (§5.1–5.4).**

- `dataset.py`:
  - `Product`, a frozen dataclass, replaces the product dictionaries.
  - `DataSet` replaces the `load_setup` tuple. Its methods are `from_file`,
    `products`, `select`, `science_images`, `table_values`, `raw_files`,
    `exptimes`, `pixel_scale`, `pointing`, `mjd`, `night`, `mosaic_labels`
    and `outdir`.
  - `check_single_setup`, `check_comb_groups` and `check_mosaic_groups`
    moved here, with the product-naming helpers, unchanged.
- `image.py`: `ReducedImage` has:
  - the header (read when first needed), `read`, `check_setup` and the
    table values;
  - the WCS and pixel scale, the solution state and the hashes;
  - the source table, through `photometry(require=True)`, which checks that
    it is current;
  - `zero_point`.
- `psf.py`: `PSFModel`, with `from_fit`, `from_meta`, `circular`, `fwhm`,
  `fwhm_major`, `fwhm_minor`, `beta`, `to_photutils`, `render` and
  `to_meta`.

**`pipnick/image/io.py` (§4.6).**

- `read_reduced` and `invalid_error` replace the readers of the photometry,
  the astrometry, the mosaics, `utils/rgb.py` and `starfind.as_masked_array`.
- The reduction still reads its calibration frames with `CCDData`, because
  `ccdproc` needs `CCDData` objects.

**Pipelines.**

- No pipeline imports another pipeline (§3.2, item 1).
- Calibration and mosaics:
  - their image dictionaries are now `ReducedImage` objects;
  - the calibration's band and coordinate flags are small helpers on the
    image;
  - the mosaic labels come from `DataSet.mosaic_labels`.
- Astrometry: `AstrometryImage` keeps its explicit constructor (the solving
  state) and holds the `ReducedImage` it came from.
- `photometry.moffat.fwhm_to_gamma` and `gamma_to_fwhm` replace the
  duplicated static methods of the fit classes, and the formula in the
  astropy-model `fwhm` properties.

**Behaviour changes** (minor):

- The calibration takes each image's filter, airmass and MJD from the
  metadata table, rather than from its source table. They are the same
  unless the table was edited and the reduction not rerun, in which case
  the table, as the authority, now wins.
- An image from another setup now stops the photometry before any image is
  measured, rather than when the loop reaches it.

**Tests.**

- New: `test_psf.py`.
- The test fixtures use `DataSet`.
- The tests of the moved functions import them from
  `pipnick.datamodel.dataset`.

**Updated:** the release notes, `docs/mosaic.rst` (the major-axis FWHM) and
the API pages.

#### Step 4: logging, verdicts, and `pipnick_status` (complete, 2026-10-28, branch `status`)

The author's revisions to the plan:

- **Logs:**
  - written by default to `<rdxdir>/log/`, with a timestamp to the second;
  - never overwritten or appended to;
  - `--log_file` still sets a different path.
- **Version:** the full pipnick version is printed first and recorded in
  every product. A step warns when it reuses products made by another
  version, encouraging (but not requiring) a rerun of all the steps.
- **Range checks:** only the bias level, read noise, FWHM and astrometric
  RMS. The sky rates and zero points of §6.2 and P6 are dropped.

Decisions made with the author:

- **Version keywords:** one per step (`PNVERS`, `PNASTVER`, `PNCALVER`).
- **Log directory** for scripts without a single reduction directory: the
  nearest output directory.
- **Failures:** the reduction and mosaics continue past a failed product, and
  a failed step exits non-zero.
- **Run records** (`runs.jsonl`, §7 P3): deferred.

Checks after the work:

- `tox -r -e 3.13-test,codestyle` passes in a fresh environment: 193 tests
  passed, 1 skipped.
- The docs build without warnings.
- All nine scripts start.
- **Real data:** the new scripts were run on a scratch copy of the 2026-08-17
  n3 (fast) reduction.
  - Each run wrote a new log, starting with the version.
  - Every product made before this step was reported as made by an earlier
    version (none recorded), with the rerun recommended, but none was
    rebuilt.
  - The reduction's progress was sequential.
  - The photometry checked every FWHM, including the 2.85″ of M11 V C096
    (inside 0.8–3.5″).
  - A dry-run astrometry gave its verdict without submitting anything.
  - `pipnick_status` reported all 19 images, and left the metadata file
    unchanged.

**Logs (`scripts/scriptbase.py`, `pkg/logger.py`).**

- `init_log(args, log_dir)` writes
  `<log_dir>/<script>_<YYYYMMDDTHHMMSS>.log`. A counter is added if the name
  exists.
- Files are opened in exclusive-create mode, and an existing `--log_file`
  is an error.
- The first lines are the version and the log path.
- Each script finds its reduction directory before logging starts, with
  `datamodel.dataset.reduction_dir`, which never changes the table.
- `entry_point` exits with the status returned by `main`.

**Verdicts (`pkg/report.py`).**

- Each step fills a `StepReport` as it runs, with failures, problems, notes,
  the versions of the reused products, and the range checks.
- `finish` logs `<STEP> OK`, `NEEDS ATTENTION` or `FAILED`, repeating the
  problems, then the version warning (once) and the next command.
- Each `*_all` takes `report=None`. The scripts pass one and return its exit
  status.

**Versions.**

- `PNVERS` is written by the reduction and the mosaics (and the group WCS
  file), `PNASTVER` by the astrometry, and `PNCALVER` by the header
  calibration.
- `pipnick_version` is in the metadata table, the source tables, the
  calibration TOML, the `_cal.ecsv` and residual tables, and every summary.
- None of these is in a hash.

**Expected ranges.**

- `[expected]` is a new section in both camera TOMLs: `readnoise_tolerance`
  0.25, `fwhm` 0.8–3.5″, `astrometric_rms` ≤ 0.5″. The bias level reuses
  `[quality] bias_level`.
- `NickelCamera.expected_range` returns each range; `StepReport.check_range`
  applies it.

**Reduction.**

- A run header gives the frame counts and products, then `[n/N]` progress in
  processing order.
- One line per frame, naming it; the details are at DEBUG.
- A failed product is recorded in the summary, and the reduction continues.
- `PNBIASLV`, `PNRDNOIS`, `PNLEVEL`, `PNSKY` and `PNCRFRAC` are saved to the
  headers, so the measurements can be reported for products that are up to
  date.
- `reduce/summary.ecsv` gives each product's status and measurements.
- The verdict collects a missing bias, the filters without flats, frames
  combined with offsets, outputs left from an earlier reduction, and the
  bias-level and read-noise checks.

**Other steps.**

- **Astrometry:** detection progress, the RMS checks, and the unsolved images
  in the verdict.
- **Photometry:**
  - progress, and the FWHM checks;
  - any exception in an image is recorded as a failure, and the photometry
    continues;
  - the jargon from `starfind` and `psf_photometry` is at DEBUG.
- **Calibration:** the failed bands and the extinction fallbacks need
  attention.
- **Mosaics:**
  - progress, and a failed mosaic or WCS group is recorded while the others
    are constructed;
  - the Gaussian-fallback warning names the image.
- **`pipnick_metadata`:** non-zero exit status on invalid arguments, a
  next-step hint, and a note when the checks that need the pixel data were
  skipped.

**`pipnick_status` (`pipelines/status.py`, `scripts/status.py`).**

- One row per science image, with the states of the reduction, astrometry,
  photometry, calibration and mosaic (`ok | stale | missing | failed | n/a`),
  the RMS, number of sources, FWHM and zero point, and an `attention`
  column.
- It ends with the next command.
- The supporting data-model additions:
  - `DataSet.from_file(read_only=True)`;
  - `DataSet.calibrations`, the bias and flat the reduction would use;
  - `ReducedImage.reduction_status`, the reduction's currency check.

**Tests.**

- New: `test_report.py`, `test_scriptbase.py`, `test_status.py`.
- New cases in `test_reduction.py` (the summary, measurements and verdict;
  a version mismatch; continuing past a failure; the reduction status) and
  `test_cameras.py` (the expected ranges).
- The version is asserted in the tests of every step.

**Docs:** a new `docs/status.rst` (logs, versions, verdicts, expected ranges,
`pipnick_status`), the release notes, and the API pages.

**Not done:**

- The run records (`runs.jsonl`), deferred.
- The rewrite of `docs/basics.rst` (P7), deferred with the other usage
  documents.
- The status vocabularies of the steps' own summaries are unchanged;
  `pipnick_status` maps them to its common vocabulary.

#### Step 5, Phase A: restructuring (complete, 2026-10-28, branch `qa`)

Decisions made with the author before starting:

- Step 5 is done in three phases on one branch, with a check-in after each:
  A, the restructuring (no change in the outputs); B, QA by default, the
  reduction QA page and the P7 fixes; C, `pipnick_report` and
  `pipnick_inspect`.
- `metadata` becomes a package, and every caller imports from its modules;
  there are no re-exports.
- (For Phase B) the removed sources go in a separate file, and the report
  uses relative links.

Checks after the work:

- `tox -r -e 3.13-test,codestyle` passes in a fresh environment: 198 tests
  passed, 1 skipped.
- The docs build without warnings.
- All nine scripts start.
- No pipeline imports `matplotlib` or another pipeline.
- **End to end:** on scratch copies of the 2026-08-17 n3 (fast) reduction,
  with the M11 images labelled as four mosaics in one WCS group, `main` and
  this branch ran the reduction, photometry (remeasured), calibration (zero
  point, applied to the headers and catalogues, with QA) and mosaics (with
  QA).
  - All reduced images (data, masks, errors and headers) are identical.
  - The source catalogues differ, but by as much as two runs of `main`
    differ from each other (9–1,800 lines per catalogue): **the photometry
    is not deterministic from run to run.** This predates the branch; it
    should be explained (and, if possible, removed) in a later step.
  - The differences in the calibration (residuals, `_cal.ecsv`) and in the
    I, R and V mosaic data (their flux scales) follow from the photometry. The other
    differences are the calibration's creation time and hash, the hashes
    that include the copies' absolute paths, and the mosaic QA figures,
    which now hide masked pixels.
  - Run without `--overwrite` on its own outputs, and on the reduction and
    photometry made by `main`, this branch reports every product as up to
    date and rebuilds nothing.

**New modules.**

- `pkg/summary.py`: `Summary`, the step summary tables (declared columns,
  types and formats; masked missing values; the version in the metadata).
  Used by every step.
- `qa/`: the figures of every step, moved out of the pipelines and the
  algorithm modules: `display` (`show_image`, `save_figure`), `sources`,
  `psf`, `calibration` and `mosaic`. Only `qa/` imports `matplotlib`.
- `datamodel/catalog.py`: `SourceCatalog`, the `_phot.ecsv` tables (column
  access, the PSF model, the celestial coordinates, the usable sources,
  instrumental magnitudes, and the residual-detection neighbours). It
  replaces `pipelines.photometry.celestial_coordinates` and
  `calibrate._instrumental` and `neighbors`.
- `datamodel/calibration.py`: `Calibration`, the night's calibration TOML
  (read, write, edits, bands, summary rows), with `toml_value`,
  `calibration_hash` and `extinction_source`.

**Module splits** (definitions moved unchanged; private names used in more
than one module were made public).

- `astrometry/qa.py` → `astrometry/header.py` (writing and removing the
  solution), `cache.py` (the astrometry.net cache and its invalidation) and
  `solve.py` (the pointing groups and the solving session).
- `photometry/fit.py` → `photometry/standards.py` (matching the standards),
  `nightfit.py` (the night's fit and its checks) and `apply.py` (applying
  the calibration to the headers and catalogues). `landolt_band` was added
  to `photometry/landolt.py`.
- `reduction/processing.py` (the frame processing and combination) and
  `reduction/combiner.py` (moved from `utils/`).
- `mosaic/prepare.py`: preparing the images and measuring their offsets.
- `metadata.py` → the `metadata` package: `table`, `setups`, `groups`,
  `names` and `build`.

**Pipelines** (lines after the work): astrometry 524, calibrate 541, mosaic
803, photometry 472, reduction 814.

**Tests:** new `test_summary.py` and `test_catalog.py`; the other tests use
the new modules.

**Updated:** the release notes and the API pages.

#### Step 5, Phase B: QA by default and small output fixes (complete, 2026-10-28, branch `qa`)

Decision made with the author:

- **Missing values in FITS tables:** the mosaic `INPUTS` columns are masked,
  so their missing values are written as NaN, the FITS null for floats, and
  masked again when read. This is the one exception to the no-NaN rule,
  limited to masked float columns of FITS binary tables.

Checks after the work:

- `tox -r -e 3.13-test,codestyle` passes in a fresh environment: 202 tests
  passed, 1 skipped.
- The docs build without warnings.
- **Real data:** on a scratch copy of the 2026-08-17 n3 (fast) reduction
  (with three existing figures deleted, and the M11 images labelled as four
  mosaics in one WCS group), every step was run twice with the defaults; the
  astrometry only as a dry run, so nothing could be submitted.
  - First run: the 32 missing figures were drawn (24 reduction products, 1
    astrometry, 2 photometry, 1 calibration, 4 mosaics), with nothing
    rebuilt. Redrawing a missing photometry figure also redraws its pair.
  - Second run: no figure was redrawn.
  - Remeasuring one image redrew its two figures and wrote its
    `_removed.ecsv` (263 sources: 96 moved, 119 non-positive, 39 outside,
    9 failed; no NaN).
  - Drawing the missing figures took about 1 s per reduction product and
    about 5 s per photometry image (restacking the PSF and subtracting the
    background).

**QA by default.**

- Every step has `qa=True` by default, and every script `--no_qa` in place of
  `--qa`.
- A figure is redrawn when its product is rebuilt, and drawn for any product
  that lacks one. For the astrometry, "rebuilt" means solved in this run, so
  a dry run only draws the missing figures.
- Figures without the product's inputs in memory are drawn from the files:
  the photometry restacks the PSF stars from the image with the saved
  settings (`build_moffat_psf_model`), the calibration reads its residual
  table, and the mosaics read the mosaic file.
- A figure that fails needs attention, but does not stop the step.
- `pkg.paths.qa_figure(outdir, name)` names every figure
  (`<step>/qa/<name>.png`), for the report in Phase C.

**Reduction QA.** `qa/reduction.py` (`plot_product`): each product, with
the histogram of its pixel values and a title with its measurements, in
`reduce/qa/<product>.png`.

**P7 fixes.**

- **Removed sources:** `vet_sources` also returns the removed sources, which
  are written to `photometry/<image>_removed.ecsv` (`pkg.paths.removed_table`,
  `ReducedImage.removed_file`): the detected and fitted positions, the
  fitted flux (masked if the fit failed), `iter_detected`, and the reason.
  Existing source tables gain the file when they are next measured.
- **Calibration summary:** new `zp` and `zp_err` columns (the zero point
  above the atmosphere, for color 0), and `file` renamed `action`.
- **Mosaic `INPUTS`:** `airmass`, `fwhm_in`, `fwhm_out` and `match_ee` are
  masked instead of -1; the QA figure uses the masks.

**Tests:** new `test_reduce_qa`, `test_photometry_qa`, `test_vet_sources`
and `test_qa_figure`; the QA, removed-source, summary and `INPUTS` checks
added to the astrometry, photometry, calibration, mosaic and path tests.

- With QA on by default, the suite took 193 s, against 78 s with QA off.
  Every pipeline call in the tests now passes `qa=False`, except in the six
  tests that check the figures (`test_reduce_qa`, `test_qa_and_detection`,
  `test_photometry`, `test_photometry_qa`, `test_calibrate_all`,
  `test_mosaic`), which also check that each figure is a readable PNG.
  The suite takes 105 s.
- The scripts' `--no_qa` sets `qa` (`dest='qa'`, `action='store_false'`).

**Docs:** a "QA figures" section in `docs/status.rst`; the `--no_qa`
options, the removed-source table, the calibration summary columns and the
`INPUTS` nulls in the step pages; the release notes and the API pages.

**Not done:** the rewrite of `docs/basics.rst` (P7) stays deferred with the
other usage documents.

#### Step 5, Phase C: `pipnick_report` and `pipnick_inspect` (complete, 2026-10-28, branch `qa`)

Decisions made with the author:

- **Shared code:** a new `pipnick/report/` package, used by the status,
  report and inspect pipelines, so that no pipeline imports another.
- **The report is written only by `pipnick_report`**, never by the steps.
- **`pipnick_inspect --display`** opens an interactive window; `--ofile`
  writes the figure instead.
- **Mosaics** are found in `<rdxdir>/mosaic`, or in `--mosaic_dir`.

Checks after the work:

- The tests pass (208 passed, 1 skipped), and the code-style check.
- The docs build without warnings.
- All eleven scripts start.
- **Real data:** on a scratch copy of the 2026-08-17 n3 (fast) reduction,
  with every step run (the astrometry only as a dry run) and the M11 images
  labelled as four mosaics:
  - `pipnick_report` wrote 25 pages in 4 s, and all 608 links resolve.
  - `pipnick_inspect` traced the standard 111 1965 in SA 111 1965 V C021 by
    `--id` and by sexagesimal coordinates (used, residual -0.023 mag, with
    its magnitudes in the seven other calibrated images), and a removed
    (`moved`) source by its coordinates.  `--ofile` wrote the display.

**`pipnick/report/`.**

- `status.py`: the status logic, moved from `pipelines/status.py` (which now
  only logs it, `status_all`).  The per-step state functions are public, and
  the mosaics can be in another directory (`mosaic_file`, `mosaic_dir`).
- `records.py`: `Record` (sections of label-value pairs and the QA figures),
  `image_record` (the provenance of an image, step by step), `trace_source`
  (a source through the source, removed-source, calibrated and standard-star
  tables), and `format_record`.
- `html.py`: the rendering of tables, records, figures and pages, with
  relative links and the states colored.

**Pipelines and scripts.**

- `pipelines/report.py` (`write_report`) and `scripts/report.py`
  (`pipnick_report`): `<rdxdir>/report/` with `index.html`, a page per
  step, and a page per image.
- `pipelines/inspect.py` (`inspect_image`, `parse_coordinates`),
  `qa/inspect.py` (`plot_inspection`) and `scripts/inspect.py`
  (`pipnick_inspect`).
- `pipnick_status --mosaic_dir`.

**Tests:** new `test_report_inspect.py` (the report's links and pages, the
records, the traces by id, coordinates and removed source, the display to a
file, the calibration page and its fallback banner).

**Docs:** a new `docs/report.rst`; the status page, the release notes and
the API pages.

**Found while checking (not changed):**

- `pipnick_status` reports the images of a band whose calibration failed as
  `calibrated: missing`, so it suggests `pipnick_calibrate` as the next step
  even though rerunning it cannot help.
- Some sources removed as `failed` were detected (`x_init`, `y_init`) far
  outside the image by the second iteration of the PSF photometry (4 of 39
  in SA 111 1965 V C021).
