# Plan: image mosaics (`mosaic` branch)

> Approved 2026-10-23.  Implementation waits for the user's go-ahead.
> Real-data checks use the user's reductions of 2026-08-17 n2 and n3 (FITS
> only; outputs in the scratchpad).

## Context

pipnick has no way to combine reduced images.  Users need mosaics of
slightly offset pointings of the same field (dithered NGC 7331 in V, R, I on
n2), and same-band stacks of non-dithered frames to register them (M11 in
B, V, R, I on n3).  Mosaics in one *WCS group* (e.g. the B, V and R mosaics of
one object) must share an identical output WCS and shape, so that they can be
combined into colour images.  The work uses the released `reproject` package.

Findings that shape the design:

- **n3 M11 spans two setups.** It was observed with both the 0.05 MHz and the
  1.0 MHz readout, i.e. two setups, two tables and two reduction directories.
  So a single-night mosaic must be able to draw from several tables.
- **Object names are unreliable.** R-filter frames are labelled
  'NGC 7331 I' or 'M11 V', so grouping uses the filter column and the
  pointing, never the name.
- **`reproject_and_coadd` doesn't fit our data model.** It neither propagates
  uncertainties nor uses masks: invalid pixels are NaN.  pipnick therefore
  reprojects the data, a good-pixel map and the variance itself, and combines
  them, converting `reproject`'s NaN to masks immediately via the footprint.
- **The dithers overlap by only about 50%.** The Andor field is about 5.7′,
  and NGC 7331 (about 10′ × 4′) is larger than the field.  A box-64 2D
  background would remove much of the galaxy.
- **The PSF model is already stored.** The photometry tables hold the Moffat
  parameters (`meta['psf']['par']`, `fittype`, `fwhm_pix`).
- **The new columns won't trigger reprocessing.** `PNTABH` covers only the
  `table_values` quantities.

## Decisions (user)

- `mosaic` and `wcs_group` are **string labels**.  Empty means not mosaicked.
- **Not defined by default.** `pipnick_metadata` adds the columns empty.
  Automatic assignment is opt-in with a command-line option.
- **PSF matching per mosaic by default,** with an option to match the whole
  WCS group.
- **Mosaics may combine several tables (setups) in this version.**

## Design

### 1. Metadata (`pipnick/metadata.py`, `pipnick/scripts/metadata.py`)

- **The columns.** Per-setup (split) tables get the masked string columns
  `mosaic` and `wcs_group`, placed after `comb_id`.  Older tables without
  them read as empty.  Appended rows get empty values, like `comb_id`'s
  `assign=is_new` handling.
- **`--assign_mosaics`** (opt-in), for new rows only:
  - Science frames are clustered by pointing (table RA/Dec within
    `--mosaic_radius`, default half the field of view).
  - `wcs_group = <objbase>_<n>`, where `objbase` is the most common object
    name in the cluster, stripped of a trailing filter token and made safe.
    `mosaic = <wcs_group>_<filter>` (primed filters written with `p`).
  - Clusters of one image get no mosaic.
  - The labels are deterministic, so the two n3 setup tables get the same
    labels.
- **Validation** (`check_mosaic_groups`, in the style of `check_comb_groups`):
  - all frames of one `comb_id` share their labels;
  - a mosaic belongs to exactly one `wcs_group`;
  - a mosaic with mixed filters gives a **warning** (and records it).

### 2. Pipeline (`pipnick/pipelines/mosaic.py`, `mosaic_all(tables, outdir=None, ...)`)

1. **Collect the images.**
   - For each table, call `load_setup` and `plan_products`, keeping the
     science products with a `mosaic` label.
   - Generalize `calibrate._collect_images` to a shared helper, passing the
     tag of each table.
   - Require a pipnick WCS (`PNASTMTH`) and a current photometry table
     (`image_hash`, `wcs_hash`).
   - Group by `wcs_group`, then by `mosaic`, across the tables.
   - Validate across the tables:
     - the mosaic-to-group mapping is unique;
     - the images overlap, using their WCS footprints, not the header RA/Dec;
     - the images come from one night (multi-night is an error for now; see
       section 4).
   - Log the setups and nights of every mosaic.
2. **Prepare each image.**
   - Convert to e⁻/s: divide the data by the exposure time, compute
     var = UNCERT²/t², and set the good-pixel map g = ~MASK.
   - Subtract the background (`--background constant|2d|none`, default
     `constant`): a sigma-clipped median with the sources masked; `2d` uses a
     large box with `subtract_background`.
3. **Match the PSFs** (`--psf_match mosaic|group|none`, default `mosaic`).
   - **Target:** a circular Moffat with FWHM = 1.05 × the worst FWHM (in
     arcsec, using the WCS pixel scale) and the smallest β.
   - **Kernel:** `photutils.psf.matching.create_matching_kernel` between the
     Moffat of each image (from `meta['psf']`) and the target, with a
     `SplitCosineBellWindow`, normalized.  The kernel is validated (residual,
     negative fraction); if it fails, the Gaussian difference kernel is used,
     with a warning.  `--kernel moffat|gaussian` selects it.
   - **Convolution:** normalized FFT convolution, conv(d·g)/conv(g).  Pixels
     where conv(g) is below a threshold are masked.  The variance is
     convolved with K (a *variance density*, correct for apertures larger
     than the noise correlation; recorded as `PNUNCTYP='density'`).
   - **QA:** the PSF is refit on each matched image (`build_moffat_psf_model`)
     and the achieved FWHM recorded.
   - An image without a PSF model (`fwhm_pix = -1`) is an error, unless
     matching is off.
4. **Scale the fluxes** (`--scale sources|calib|none`, default `sources`).
   - Cross-match the sources one-to-one (mutual nearest neighbours) by sky
     position.  Only unflagged, unsaturated sources above an S/N cut are used.
   - With PSF matching on, measure aperture fluxes on the matched images
     (`aperture_analysis`, radius 1.5 × the target FWHM).  With matching off,
     use `flux_psf`/t.
   - Solve **globally** for the log scales from all overlapping pairs:
     weighted least squares with a 0.01 mag floor and sigma clipping, anchored
     to a reference image (lowest airmass by default; `--reference`).
   - Check that the overlap graph is connected.  Images with fewer than
     `--min_common` stars (default 3), or in an unconnected part of the graph,
     get a warning and are not scaled.
   - `calib` uses `PNCALZP` when every image has it.  If the reference has
     `PNCALZP`, the mosaic gets the corresponding `PHOTZP`.
   - Optional: additive background offsets solved in the overlaps after the
     scales (`--match_background`).
5. **Output WCS (one per WCS group).**
   - `find_optimal_celestial_wcs` over every image in the group (all bands),
     with `frame='icrs'`, north up, and the finest input pixel scale (or
     `--pixel_scale`).
   - Every mosaic in the group uses this identical WCS and shape, saved as
     `mosaic/<wcs_group>_wcs.fits` (a header only).
6. **Reproject and combine.**
   - Reproject d·g, g and var·g with the same routine (`--reproject adaptive`
     with flux conservation, the default; or `interp` or `exact`).  Correct for
     the pixel-area ratio r: flux ×r, variance ×r.
   - Mask output pixels with footprint 0, or g below 0.99, immediately.
   - Combine with **one weight per image**, 1/(s² × background variance),
     times the reprojected g, rather than per-pixel weights (those would bias
     bright sources).
   - With 3 or more contributors, a sigma-clipped weighted mean
     (`--rej`, `--no_reject`).
   - The combination is done in row chunks to limit memory.
7. **Write the products** (`mosaic/<mosaic>.fits`, by default in
   `rdxdir/mosaic` for one table; `--outdir` is required for several).
   - **Extensions:** primary (e⁻/s, `BUNIT='electron/s'`), UNCERT, MASK,
     EXPMAP (Σ exposure time) and COVERAGE (the number of images).
   - **Header:**
     - the WCS;
     - the filter, the target FWHM, the kernel type and the reference image;
     - the inputs, with their setups, nights, image/WCS/table hashes, scales
       (with errors) and achieved FWHMs;
     - the settings hash and the group-membership hash;
     - warnings, such as mixed filters or unscaled images.
   - **Summary:** `mosaic/summary.ecsv`, also logged.
   - **QA** (`--qa`): the mosaic, its coverage, the scale fits, and the PSF
     matching.
   - **Reruns:** a mosaic is rebuilt only if its inputs, the settings or the
     group membership change.  A change in membership rebuilds the whole
     group, since the WCS changes.
8. **Script** (`pipnick/scripts/mosaic.py`, entry point `pipnick_mosaic`):
   - `-t` takes one or more tables; `--outdir`, `--mosaics` (select labels),
     `--overwrite`, and the options above;
   - optional `--config <toml>` with the tables, output directory and
     settings.  Membership stays in the tables.

### 3. Dependencies, docs and release notes

- `reproject` (a released version, e.g. `>=0.14`) becomes a core
  dependency; the installed version is 0.21.
- A new `docs/mosaic.rst`.  It covers the workflow (reduce → astrometry →
  photometry → [calibrate] → mosaic), the columns, the PSF matching, the
  variance-density convention, the scaling, and the multi-table use.
- API pages and release notes.

### 4. Towards multi-night mosaics (infrastructure built now)

- **Labels are global.** The same label in different tables joins, so a
  multi-night mosaic is the same mechanism as the multi-setup one.
- **Opt-in.** For now, mosaics spanning nights are an error; a later
  `--allow_multi_night` flag (or the config) enables them.  Opt-in is needed
  because default labels like `M11_V` could otherwise merge separate nights
  silently.
- **Outputs and provenance are independent of any one reduction.** Products
  can go to `--outdir`, and each input records its table, setup, night and
  hashes.
- **Nothing else assumes one night.** The scaling (global and relative)
  already handles transparency differences between nights.  The per-night
  photometric calibrations are used through `--scale calib`.
- **The config file** becomes the record of a multi-night project: the list
  of tables, the output directory and the settings.

## Reused code

- `pipnick/pipelines/reduction.py`: `load_setup`, `plan_products`.
- `pipnick/pipelines/photometry.py`: `_image_hash`, `_wcs_hash`, `_select`.
- `pipnick/pipelines/calibrate.py`: `_collect_images`, moved to a shared
  helper; `_cross_match`, adapted to one-to-one matching.
- `pipnick/photometry/psf_photometry.py`: `subtract_background`,
  `build_moffat_psf_model`.
- `pipnick/photometry/aperture_photometry.py`: `aperture_analysis`.
- `pipnick/photometry/moffat.py`: the Moffat models.
- `pipnick/astrometry/fit.py`: `wcs_from_header`.
- The `pipnick/astrometry/qa.py` figure pattern.

## Verification

- **Tests** (`pipnick/tests/test_mosaic.py`, synthetic, no network): small
  TAN-WCS images and photometry tables written directly.
  - A Moffat kernel equalizes the encircled energy to better than 1%; the
    Gaussian fallback does worse.
  - Normalized convolution keeps a constant image constant.
  - A Monte Carlo of aperture-sum variances agrees with UNCERT to within
    about 10%.
  - Injected scales are recovered to better than 0.5% over a chain of dithers
    with no star common to the whole chain.  Too few common stars gives a
    warning and no scaling.
  - Flux is conserved, the outputs are finite (`np.isfinite`), and the masks
    are right.
  - An injected cosmic ray is rejected when there are 3 or more images.
  - Metadata:
    - the columns are empty by default;
    - the opt-in clustering;
    - the `comb_id` consistency;
    - the mixed-filter warning;
    - a two-table join with the same labels;
    - the identical WCS and shape within a group.
  - Reruns: unchanged inputs are a no-op; a change in membership rebuilds
    the group.
- **Checks:** `tox -e 3.13-test,codestyle`, and the docs build without warnings.
- **Real data** (scratchpad):
  - n2 NGC 7331, V/R/I in one WCS group: check the seams, the scales and the
    matched FWHM.
  - n3 M11, per band across both setup tables: check the registration, and
    that the shared WCS gives aligned B, V and R mosaics (a quick 3-colour
    PNG in the scratchpad).

## Implementation progress (2026-10-26)

Implemented:

- the metadata columns and `--assign_mosaics`;
- `pipnick/mosaic/` (`psfmatch`, `scaling`, `coadd`);
- `pipnick/pipelines/mosaic.py` and `pipnick_mosaic`;
- `reproject` as a core dependency;
- `pipnick/tests/test_mosaic.py` (12 tests);
- `docs/mosaic.rst`, the API pages and the release notes.

`tox -e 3.13-test,codestyle` passes, and the docs build without warnings.
`_safe_name` and `_filter_name` moved from `reduction` to `metadata` (reduction
imports them), so that the labels can use them.

Deviations from the plan, found on the real data:

- **Target power index:** the target PSF uses the *median* Moffat β of the
  images, not the smallest.  The β of the photometry PSF models is noisy
  (1.1–1.6 on n2; at the fit bound of 6 on M11), and the smallest selects the
  noisiest measurement.
- **Kernel taper:** the Moffat kernel is tapered to a circle, which removes
  the square halos around bright stars.  It costs up to about 2% in encircled
  energy for very low-β targets.  The accuracy of each image's kernel is
  recorded in the `INPUTS` table (`match_ee`).
- **Kernel acceptance:** a kernel is rejected only if it amplifies the noise
  (Σ K² ≥ 1), not for its negative values.  Near-delta kernels ring, but are
  accurate.
- **QA refit:** the matched FWHM is measured by a refit with β fixed at the
  target value.  Free refits are degenerate.
- **Local background for scaling:** the scaling apertures use a local
  background annulus.  The global background was biased by the stars' wings.
- **Rejection:** the rejection compares each image with the mean of the
  *others* and rejects at most one image per pixel per iteration, with a
  fractional term (`--rej_frac`, default 0.1).  Otherwise, the cores of bright
  sources were rejected.
- **New QA:** an `NREJECT` extension, and a warning for seeing outliers.
- **Background matching:** `--match_background` is implemented (off by
  default).  On n2, it removes the seams between dithers.
- **Memory:** images are reprojected into their bounding boxes rather than
  combined in row chunks.
- **`--reference`:** applies only to the mosaic that contains that image.

Real data (scratchpad only):

- **n2 NGC 7331:** V, R and I in one WCS group; rebuilding is a no-op.
- **n3 M11:** combined across both setups.  The matched FWHM is on target
  (3.0″), and the scale errors are about 0.002 mag.  The 3-colour R/V/B image
  with `--psf_match group` is aligned.

Open items:

- **n3 append fails:** `pipnick_metadata --append` on n3 raises an error.
  The raw directory has 98 frames not in the user's tables, and one of the
  setups they form mixes image widths (1044 and 2088).  This is a
  pre-existing per-setup split issue, not specific to the mosaics.
- **Photometry fit failure:** the photometry Moffat fit fails on broad-PSF
  images ("initial guess outside bounds"), so the QA FWHM is missing for
  M11_V_C096.  This is recorded with the photometry follow-up.
- **Residual rejections:** a few regions still differ by more than 20% between
  images (near-saturated cores, or a poorly matched image).
- **Matched FWHM on n2:** the measured FWHM is above the target in the
  NGC 7331 field, which has few clean stars.  On M11 it is on target.
- **Multi-night mosaics:** not enabled yet.

### Colour images and single-image mosaics (2026-10-26)

- Added `pipnick_rgb` (`pipnick/scripts/rgb.py`, `pipnick/utils/rgb.py`).
  - It makes a three-colour image from three images with the same WCS and
    shape: Lupton et al. (2004) by default, or a linear or log stretch.
    Each channel has its background subtracted and is normalized by a
    percentile.
  - Images without a shared WCS raise an error that recommends
    `pipnick_mosaic` with a shared `wcs_group`.
- Mosaics with a single image work.  With `--psf_match mosaic` they are not
  convolved, only reprojected.
- Tests: `test_rgb.py` (2) and `test_single_image_mosaics`.  Docs and release
  notes are updated.
- Real data: the single-image R/V/B mosaics of M11 (n3, fast readout) register
  well with V image C099.  With V image C096, V is offset by about 3 pixels.
  That frame looks bad: its astrometric rms is 0.18″ (others 0.03″), from
  only 43 matches, and its seeing is 2.86″.
