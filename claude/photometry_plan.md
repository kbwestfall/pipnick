# Assessment and plan: update the photometry pipeline (`photometry` branch)

> Design record for the `photometry` branch. Status: implemented and verified
> on synthetic and real data; awaiting review. Progress is recorded at the end
> of this file.

## Context

`pipnick_photometry` measures the sources in the reduced science images. It
predates the metadata, reduce and astrometry branches and does not run.

The reduction now writes `<object>_<filter>_C<comb_id>_<date>.fits` images
with `MASK` and `UNCERT` extensions, `PN*` provenance, and table-authoritative
object, filter, exposure time and coordinates. The astrometry adds solved WCS
headers (`PNAST*`).

**User requirements for this branch:**
- Photometry must run without an astrometric solution.
- Pixel coordinates are always provided.
- Celestial coordinates are provided if the image has a WCS solution, and are
  masked otherwise.

## Assessment of the current code

### `pipnick/pipelines/photometry.py` (`photometry_all`)

- **Mid-refactor, not runnable.** It stops at `embed(); exit()` (lines
  131-132) before consolidating, aperture photometry, or writing. It writes to
  an undefined `output_file`.
- **Finds images by the old naming:**
  - `metadata.default_reduced_file_name` (`<raw>_rdx.fits`);
  - `build_camera_datatable(..., ext='_rdx.fits')`.

  Neither matches the reduction's outputs.
- **Error handling:**
  - `return ValueError(...)` (line 102) returns the exception instead of
    raising it.
  - The `phot_table` output is set up but never written.

### `pipnick/scripts/photometry_script.py`

- Mutually exclusive `-r`/`-t` options, and `-p` (`phot_table`).
- Large blocks of commented-out retired options remain.
- The file has not been renamed to `photometry.py` (planned).

### `pipnick/photometry/` modules

| Module | Role | State |
|---|---|---|
| `starfind.py` | Bright-star selection for the PSF model (`detect_sources`, `kron_stats`, `extract_stamps`, `plot_stamps`) | **Broken on photutils 3:** uses `sources['xcentroid']`, which is now `x_centroid`. Otherwise reusable. |
| `moffat.py` | Moffat fitters (`FitMoffat2D`, `FitEllipticalMoffat2D`) and the astropy model `EllipticalMoffat2DModel` | Self-contained; looks sound. |
| `psf_photometry.py` | `build_moffat_psf_model`, `psf_analysis`, `consolidate_groups` | **Broken on photutils 3**; see below. |
| `aperture_photometry.py` | `aperture_analysis`, `format_table` | **Broken:** the signature takes `x, y`, but the body uses an undefined `phot_data`. Still uses `Fits_Simple`. |
| `fit.py` | Legacy stamp fitting (`fit_psf_stack/single/generic`) and plotting (`psf_plot`, `plot_sources`) | `fit_psf_*` read stamp files that nothing writes any more (`*.rdx.fits` `STAMPS_*`, `*.src.db`), so they are dead. Only `psf_plot` and `plot_sources` are imported. |
| `calibrate.py` | Photometric calibration (zero point, airmass and color terms, iterative rejection) | Standalone; not wired into any pipeline. Out of scope for this branch. |

**photutils 3 incompatibilities in `psf_analysis`** (checked against the
installed photutils 3.1.1):
- the keyword `localbkg_estimator` is now `local_bkg_estimator` (a
  `TypeError`);
- `IRAFStarFinder(minsep_fwhm=...)` is deprecated (now `min_separation`);
- the output column `npixfit` is now `n_pixels_fit` (used by
  `consolidate_groups` and `format_table`), and there is a new column
  `reduced_chi2`.

The mapping of the Moffat parameters to `MoffatPSF` is correct: our `gamma`
(core width) becomes `alpha`, and our `alpha` (power index) becomes `beta`.

### Related code

- **`pipnick/pipelines/final_calib.py`** adds RA/Dec to photometry catalogs
  from the old `astrometric/` files. It is superseded by the requirement
  above, which puts the coordinates in the photometry tables.
- **`pipnick/utils/fits_class.py`** (`Fits_Simple`), `display_fits.py` and
  `graphs.py` are the legacy image wrapper and plotting helpers.
  `graphs.plot_sources` duplicates `fit.plot_sources`.
- **`metadata.default_reduced_file_name`** is used only by the photometry
  pipeline.
- `from IPython import embed` is left in `photometry.py`, `psf_photometry.py`
  and others.
- `docs/photometry.rst` is stale. There are no photometry tests.

## Summary of the current algorithm

For each reduced image:

1. **Build the PSF model** (`build_moffat_psf_model`).
   - Detect sources with `DAOStarFinder` (`starfind.detect_sources`).
   - Measure each source's Kron radius and flux, with a local background
     (`kron_stats`).
   - Keep isolated sources (no neighbour within 5 FWHM) that are bright
     (Kron flux above a minimum).
   - Extract background-subtracted stamps (window 10 times the mean Kron
     radius), and keep those with less than 1% masked pixels.
   - Average the stamps, and fit a circular (or elliptical) Moffat to the
     stack. Its FWHM sets the scales below.
2. **PSF photometry** (`psf_analysis`) with photutils
   `IterativePSFPhotometry`:
   - finder: `IRAFStarFinder` at `thresh` times the MAD background RMS;
   - grouper: `SourceGrouper`, minimum separation 2 FWHM;
   - local background: `MMMBackground` in an annulus from 3 to 6 FWHM;
   - PSF model: `MoffatPSF` with its shape fixed from the stack fit (flux and
     position free);
   - fit window about 2 FWHM; aperture radius 1.8 FWHM for the initial
     fluxes;
   - repeated detection and fitting on the residual image (`mode` `'all'`
     or `'new'`).

   The header `AIRMASS` is added as a column.
3. **Optional consolidation** (`consolidate_groups`): merge the sources fit
   together in a group into one (fluxes summed, positions flux-weighted).
4. **Aperture photometry** (`aperture_analysis`):
   - circular apertures (radius 8 pixels) at the fitted positions;
   - the local background subtracted with `ApertureStats`;
   - `flux_aper` is masked for sources with masked pixels;
   - the ratio of PSF to aperture flux is added.
5. **Write** a CSV catalog per image. This step was never reached.

Photometric calibration (`calibrate.py`) and the RA/Dec conversion
(`final_calib.py`) are separate and not connected.

## Lingering TODO comments (and similar issues)

| Location | Note | Plan |
|---|---|---|
| `fit.py:189`, `fit.py:202` | `_moffat_p0` amplitude taken from row `i` of the stamp | Moot: the dead `fit_psf_*` functions are removed (decision 3). |
| `starfind.py:228` | Consider `photutils.segmentation` for detection | Defer. Note it as an alternative in the docstring. |
| `starfind.py:288` | Clipping parameters; clip high values preferentially | Expose `sigma` as a parameter; defer the rest. |
| `psf_photometry.py:70` | Mean or median stack of stamps | Add a parameter (default mean, as now), with a test of both. |
| `psf_photometry.py:136` | Make the grouper separation a parameter | Do it (`group_sep`, in FWHM). |
| `psf_photometry.py:148` | Make the elliptical Moffat a proper PSF model like `MoffatPSF` | Defer. Keep `make_psf_model` and test that it runs. |
| `calibrate.py:311` | Only need each color (magnitude pair) to be valid | Defer; calibration is out of scope. |
| Not marked as TODO | `embed()`/`exit()`; `return ValueError`; `process_par` logs an error without raising; mutable default `preserve=[]`; `LevMarLSQFitter` (legacy); `AIRMASS` read from the header instead of the table | Fix all. |

## Development plan

0. **Airmass from the start of the processing** (user request).
   - **Current state:**
     - The metadata table already has an `airmass` column from the first
       step (`pipnick_metadata`), filled from the `AIRMASS` header keyword
       for both cameras (`camera.metadata`). The real C2 n1 and Andor n1
       tables have no missing values.
     - The airmass is *not* treated like the other table values in the
       reduction. `reduction._table_values` sets the object, filter,
       exposure time and RA/Dec from the table, but not the airmass.
     - A combined science image keeps the first frame's header airmass.
     - An edited airmass in the table is not propagated, and does not
       trigger reprocessing.
   - **Changes:**
     - **Metadata:** if a frame's header has no airmass, compute it from the
       table RA/Dec, the frame's mid-exposure time and the Lick Observatory
       location (stored in the camera configuration, so that no network site
       lookup is needed). Record it in the `note` column. The table stays the
       editable authority.
     - **Reduction:** add `airmass` to `_table_values` for science images: the
       mean airmass of the frames in a group. It is written to the camera's
       airmass keyword (and `AIRMASS`), with `ORIGAIR` if it differs, and
       included in `PNTABH`, so that edits trigger reprocessing.

       Because the hash changes, existing reductions are reprocessed once on
       their next run. Calibration frames do not need the airmass.
     - **Photometry:** take the airmass from the metadata table, i.e., the
       value now also in the reduced header, and record it in each catalog's
       metadata and as a column.
   - **Tests:**
     - an airmass computed for a header without one;
     - an edited airmass reprocesses the science image and appears in its
       header;
     - the mean airmass of a combined group.

1. **Inputs like the other pipelines.**
   - `photometry_all(table, rdxdir=None, ...)` uses `reduction.load_setup` and
     `reduction.plan_products` to find the science images, skipping missing
     files with a warning and checking `PNSETUP`.
   - It reads each image with its `MASK` and `UNCERT` extensions.
   - It takes object, filter, exposure time, airmass and MJD from the
     metadata table, which is the authority.
   - `--images` selects images, as in the astrometry.
2. **Fix the modules for photutils 3, without changing the algorithm:**
   - `starfind`: column names;
   - `psf_photometry`: `local_bkg_estimator`, `min_separation`, the column
     names, `TRFLSQFitter`; and pass `UNCERT` as `error`, so that `flux_err`
     is meaningful;
   - `aperture_photometry`: a correct signature and body; work on the masked
     data; add an aperture flux error from `UNCERT`; drop `Fits_Simple`.
3. **Scales from the camera.** The initial FWHM guess is the camera's
   `fwhm_guess` divided by the pixel scale, not a fixed 7 pixels. The
   aperture radius becomes a multiple of the measured FWHM (currently a fixed
   8 pixels).
4. **Output source table per image:** `rdxdir/photometry/<image>_phot.ecsv`.
   ECSV keeps masks and metadata; CSV does not.
   - **Columns:**
     - `id`, `group_id`, `group_size`, `iter_detected`;
     - **`x`, `y`** (0-based pixel, always present) with `x_err`, `y_err`;
     - **`ra`, `dec`** (deg, ICRS; masked if there is no solution);
     - `flux_psf`, `flux_psf_err`, `flux_aper`, `flux_aper_err` (electrons;
       masked where invalid);
     - `mag_inst`, `mag_inst_err`, where `mag_inst` = -2.5 log10(flux_psf /
       exptime);
     - `local_bkg`, `qfit`, `cfit`, `reduced_chi2`, `flags`.
   - **Metadata:** image, setup, comb_id, object, filter, exptime, airmass,
     mjd; the PSF model (type, parameters, FWHM in pixels and arcsec); the
     detection settings; the image provenance hashes (`PNRAWHSH`, `PNTABH`);
     and the solution provenance (`PNASTMTH`, `PNASTWH`), or a note that
     there is none.
5. **Celestial coordinates only from a pipnick solution.**
   - The reduced images also carry the *raw-header* WCS, which is the
     telescope's nominal one and unsolved, so a WCS counts only if the header
     has `PNASTMTH`.
   - `ra`/`dec` are computed with `astrometry.fit.wcs_from_header` and
     `to_icrs` (ICRS, as for the astrometry).
   - Otherwise they are masked, and the log notes that there is no solution.
6. **Reruns.**
   - An image is re-measured only if it changed (`PNRAWHSH` or `PNTABH`), the
     photometry settings changed (hash), or `--overwrite` is set.
   - If **only the WCS** changed (`PNASTWH` differs, or a solution appeared),
     only `ra`/`dec` are recomputed. This is cheap, so running the astrometry
     after the photometry just updates the coordinates.
7. **Summary.** `rdxdir/photometry/summary.ecsv`, also logged, gives per
   image: the number of sources, the FWHM, whether there is a WCS, the
   status, and a reason if it failed.
8. **Script.**
   - Rename to `pipnick/scripts/photometry.py`, as for the astrometry, and
     update the entry point and API docs.
   - Options: `-t` (required) and `-r`; `--images`; `--overwrite`;
     `--thresh`, `--fwhm`, `--fittype {circ,ellip}`, `--mode {all,new}`,
     `--group`, `--group_sep`, `--aper_radius`; `--qa` (a figure with the
     sources and PSF fit per image).
   - The retired, commented-out options and `-p` are removed.
9. **Retire superseded code** (decisions below):
   - `pipelines/final_calib.py`;
   - the dead `fit.py` `fit_psf_*` functions;
   - `metadata.default_reduced_file_name`;
   - `Fits_Simple` uses in the photometry modules.
10. **Tests** (`pipnick/tests/test_photometry.py`; no network):
    - Synthetic reduced images built with a known Moffat PSF at known
      positions and fluxes, with noise, `MASK` and `UNCERT`.
    - **Recovery:** the PSF FWHM, positions to within about 0.05 pixel, and
      fluxes within their errors; aperture and PSF fluxes agree for isolated
      stars.
    - **Coordinates:** `ra`/`dec` are correct with a `PNAST*` solution, and
      masked for an image with only a raw-header WCS.
    - **Reruns:** a rerun skips unchanged images; a new solution updates only
      the coordinates; changed settings or a changed image trigger
      re-measurement.
    - **Other:** `consolidate_groups`; `--images`; the summary; a regression
      test that the photutils 3 column names are used.
11. **Docs.** Rewrite `docs/photometry.rst` for the new interface (outputs,
    columns, the WCS behaviour, reruns), and add release notes.
12. **Real data** (scratchpad only, reading FITS files only, no astrometry.net
    submissions): run on the C2 2024-06-24 n1 reduction, which has solutions
    for all 22 images, and on a copy with the solutions removed.
    - Check that `ra`/`dec` are filled or masked as expected.
    - Check that the PSF FWHM is consistent with the astrometry's seeing.
    - Check that a rerun is a no-op.

## Decisions to confirm

1. **Output format and place:** a per-image ECSV in `rdxdir/photometry/`,
   plus `summary.ecsv`. This replaces per-image CSV and the `-p` table.
2. **Retire `final_calib.py`:** its job (adding RA/Dec) moves into the
   photometry tables.
3. **Remove the dead code.** Not all of `fit.py` is dead.
   - **Dead (remove):**
     - `fit_psf_stack`, `fit_psf_single` and `fit_psf_generic` (with the
       `_moffat_p0` TODOs at lines 189 and 202): they read stamp files that
       nothing writes.
     - `_moffat_p0`, which only `fit_psf_generic` uses.
     - `get_graphable_pars`, which calls `generate_stamps`. That function no
       longer exists (its import is commented out), so the call would raise a
       `NameError`.
     - The commented-out `get_source_pars`.
     - `metadata.default_reduced_file_name`.
   - **Live (keep):**
     - `psf_plot` plots the stacked PSF, the Moffat model, the residuals and
       the radial profile from a `FitMoffat2D`/`FitEllipticalMoffat2D`
       object, which is exactly what `build_moffat_psf_model` produces. It
       becomes the PSF panel of `--qa`. Two fixes:
       - remove the guard that raises for `fittype='circ'` (the body already
         handles both);
       - use matplotlib `Figure` instead of `pyplot`.
     - `get_param_list` computes FWHM, eccentricity and angle from the
       Moffat parameters, a general utility for PSF analyses. Keep it.
     - `plot_sources` and `get_apertures`: an interactive display of the
       sources (used by `psf_analysis(show=True)`). Keep them; the `--qa`
       figure is the non-interactive counterpart.
       `utils/graphs.plot_sources` duplicates them; it is left for the
       future cleanup of `utils/`.
4. **Keep `consolidate_groups` and `--group`** as an option. The default
   stays off.
5. **Instrumental magnitudes** per second (`-2.5 log10(flux / exptime)`).
   Calibration (`calibrate.py`) stays out of scope for this branch.
6. **Airmass of a combined image:** the mean of its frames' airmasses
   (alternatively, the airmass at the mean time of observation).
7. **Airmass fallback: use astropy** (user decision).
   - The airmass is `SkyCoord(ra, dec).transform_to(AltAz(obstime=t,
     location=site)).secz`, with `t` the mid-exposure time.
   - `site` is an `astropy.coordinates.EarthLocation`.
   - **`EarthLocation.of_site('Lick Observatory')` cannot be used directly:**
     it downloads astropy's site registry from astropy-data on first use, and
     there is no bundled copy. Tested here with the network blocked, it fails
     (`URLError`). It would therefore fail in the tests (the network is
     blocked) and offline.
   - So the observatory's geodetic coordinates (longitude, latitude, height)
     are stored in an `[observatory]` section of the camera configuration, and
     the location is built with `EarthLocation.from_geodetic`.
   - A test marked `remote_data`, skipped by default and never run in CI,
     checks that these values agree with `EarthLocation.of_site('Lick
     Observatory')`.
   - The transform runs with astropy's IERS automatic download disabled
     (`iers.conf.auto_download = False`) and degraded accuracy allowed, so it
     never contacts the network. The resulting error in UT1-UTC (< 1 s) is
     negligible for the airmass. Tested offline: it works without warnings
     for a 2024 date.
   - The computed airmass is checked against the header `AIRMASS` of the real
     C2 frames during development.

## Critical files

- **Modified:**
  - `pipnick/metadata.py`, `pipnick/cameras/cameras.py` and the camera TOML
    files (airmass fallback and observatory location);
  - `pipnick/pipelines/reduction.py` (airmass as a table value);
  - `pipnick/pipelines/photometry.py` (rewrite);
  - `pipnick/photometry/psf_photometry.py`, `starfind.py`,
    `aperture_photometry.py`;
  - `pipnick/scripts/photometry_script.py`, renamed to `photometry.py` with
    `git mv` (plus the entry point in `pyproject.toml` and the API docs).
- **Removed, if confirmed:** `pipnick/pipelines/final_calib.py` and parts of
  `pipnick/photometry/fit.py`.
- **Reused:**
  - `reduction.load_setup` and `plan_products`;
  - `astrometry.fit.wcs_from_header` and `to_icrs`;
  - `astrometry.qa`-style figures (matplotlib `Figure`);
  - the summary pattern from `pipelines/astrometry.py` (`_summarize`).
- **New:** `pipnick/tests/test_photometry.py`; docs and release notes.

## Verification

- `tox -e 3.13-test,codestyle` passes.
- The synthetic tests recover the injected PSF, positions and fluxes.
- Real C2 n1 data, with and without WCS solutions:
  - catalogs are written with the expected columns and masking;
  - a rerun is a no-op;
  - removing and re-adding a solution updates only the coordinates.

## Progress

- **Step 0 (airmass):**
  - `[observatory]` added to both camera TOMLs, using the values of astropy's
    `EarthLocation.of_site('Lick Observatory')`: -121.636667, 37.343333,
    1290 m.
  - `NickelCamera.location` (`EarthLocation.from_geodetic`) and
    `NickelCamera.compute_airmass` (`AltAz` `secz` at mid-exposure, with the
    IERS auto-download off).
  - `pipnick_metadata` computes a missing airmass and notes it as "airmass
    computed".
  - The computed airmass agrees with the header values to 0.004 (C2) and
    0.006 (Andor).
  - Reduction: `table_values` (now public) adds the mean airmass to science
    images, with `ORIGAIR`, and includes it in `PNTABH`.
  - A `remote_data` test (skipped by default) checks the location against
    astropy.
- **Steps 1-9 (photometry):** implemented as planned, with these changes,
  found on real or synthetic data:
  - **Background:** a 2D background is now subtracted before detection; the
    original code never subtracted one.
  - **Finder:** `DAOStarFinder` replaces `IRAFStarFinder`. IRAF thresholds
    the unsmoothed image and found about 5 times more spurious sources on
    real images (109 vs 23 at 5 sigma in NGC 3982).
  - **Separation:** the minimum separation of detected sources is now 1 FWHM
    (`--min_sep`). The original 0.1 FWHM gave several detections per star.
  - **Vetting:** `vet_sources` removes failed fits, fits outside the image,
    fits that moved more than 1 FWHM, and non-positive fluxes. On real data
    this is 15-40% per image, mostly spurious iteration-2/3 detections in
    the residuals of bright stars. The counts are recorded in the table
    metadata.
  - **Pixel errors:** a model, (counts + RN^2) / NCOMBINE, instead of the
    `UNCERT` extension (whose origin is unclear for single frames).
  - **Elliptical PSF:** `photutils.psf.make_psf_model` is deprecated in
    photutils 3.1 and its compound model has no bounding box, which
    iterative mode needs. A flux-normalized `EllipticalMoffatPSF` was added,
    which resolves the TODO at `psf_photometry.py:148`.
  - **Pre-existing bugs fixed in `moffat.py`:**
    - the `FitEllipticalMoffat2D.guess_par` array had 6 elements instead of
      8, and `gamma1` used `par[4]` instead of `par[6]`;
    - the peak was indexed as `[x, y]` in both fitters;
    - `EllipticalMoffat2DModel.fwhm` used a nonexistent `gamma`.
  - **TODOs resolved:** `starfind.py` (`sigma` parameter; segmentation noted
    in the docstring), the stack mean/median parameter, `--group_sep`, and
    the elliptical PSF model. Only `calibrate.py:311` remains (deferred).
  - **Retired:** `final_calib.py`, the dead `fit.py` functions, the
    commented-out `generate_stamps`, and `default_reduced_file_name`.
    `examples/all.py` and `docs/example.rst` still reference `final_calib`;
    they are stale old examples, left for the docs rewrite.
- **Tests:** `test_photometry.py` has 10 tests, plus the airmass tests in
  `test_metadata.py` and `test_reduction.py`. `tox -e 3.13-test,codestyle`:
  116 passed, 1 skipped.
- **Real data** (C2 2024-06-24 n1, scratchpad):
  - All 22 images measured in about 2 min 20 s; PSF FWHM 1.3-1.9 arcsec.
  - `ra`/`dec` match the saved Gaia queries to a median of 0.03-0.14 arcsec,
    for all solution methods.
  - A rerun takes 2.9 s and rewrites nothing.
  - Removing a solution masks the coordinates. Re-solving it with
    `pipnick_astrometry` (0 submissions) and rerunning the photometry only
    updates the coordinates.

## Pixel errors (`UNCERT`), traced and reworked (user request)

**Origin of the old `UNCERT`:**
- pipnick never set an uncertainty; it came from ccdproc and astropy.
- Combined products got ccdproc's `Combiner.average_combine` uncertainty:
  the std of the combined values divided by sqrt(N) (`combiner.py:576-578`).
- Single-frame science images inherited only the *calibrations'*
  uncertainty, through astropy arithmetic in `subtract_bias` and
  `flat_correct`.
- On real C2 data, the single-frame `UNCERT` was about 7.5 times smaller
  than the real noise (2.8 vs 20.8 e-). The combined values (N = 2-4) were
  very noisy, biased low, and inflated at stars (243 e-, from unregistered
  frames).

**New behaviour:**
- **`process_frame`:**
  - variance = RN^2 (1 + 1/N_os) after the overscan subtraction, where N_os
    is the number of overscan columns averaged per row by
    `ccdproc.subtract_overscan` (mean per row; a nominal value, correlated
    along the row);
  - plus the bias variance (propagated by astropy);
  - plus Poisson of the bias-subtracted counts (`photon_noise`, False for
    biases);
  - the flat's error is propagated by astropy.
- **`combine_frames`:** uses `pipnick.utils.combiner.Combiner`, which
  propagates sqrt(sum (w s sigma)^2) / sum w over the frames kept after
  rejection.
- **Upstreaming:** the combiner is self-contained, with only numpy, astropy
  and ccdproc imports, and its tests (`tests/test_combiner.py`) use only
  that module. Both are ready to propose to ccdproc, where `Combiner` would
  keep the uncertainties and `average_combine` would gain
  `propagate_uncertainty`.
- **Photometry:** uses `UNCERT` directly (the noise-model function is
  removed). An image without `UNCERT` gets a warning and unweighted fits.
- **No compatibility code** for older reductions (user decision). The old
  scratchpad reductions were deleted and the night reduced again from
  scratch.

**Tests:**
- `test_combiner.py` (5 tests);
- `test_reduction.py::test_uncertainty`, which matches a hand calculation
  for the bias and a combined science image to 0.1%;
- `test_photometry.py::test_pixel_errors`, which checks that doubling
  `UNCERT` doubles the flux errors and that a missing `UNCERT` is reported.
- `tox`: 122 passed, 1 skipped.

**Real data** (C2 n1, fresh reduction `rdx_v2`):
- `UNCERT` vs the measured sky noise: median ratio 0.98 over 22 images
  (I band 0.86-0.94, likely fringing).
- Repeated exposures: MAD-based scatter of the normalized flux differences
  0.8-1.5. The errors include only pixel noise, not PSF-model or flat-field
  errors.

**Astrometry of the fresh reduction:** one approved astrometry.net submission
(ID 16336491) for `109_231_I_C020`.
- All 22 images solved (1 astrometry.net, 15 bootstrap, 6 Gaia).
- The solution is cached exactly as returned (`EQUINOX` only, no `PN*`
  keywords), with the index `astrometry_net.ecsv`.
- Rerunning the photometry updated only the coordinates of all 22 catalogs.
  They match Gaia to median offsets of 0.02-0.13 arcsec.

**Follow-up (user request): reuse astropy's propagation.**
- For unweighted combinations, `Combiner._propagated_average_uncertainty` now
  stacks the scaled, rejection-masked data and uncertainties in one 3D
  `CCDData` and takes the uncertainty of `.mean(axis=0)`
  (`NDUncertainty._propagate_collapse`).
- In a 4 x 1024 x 1024 test with rejections, this gave results identical to
  the explicit formula (maximum difference 0) and was faster (0.13 s).
- astropy has no weighted mean, so the explicit formula remains for weighted
  combinations.
- New test: `test_propagation_matches_formula`.
