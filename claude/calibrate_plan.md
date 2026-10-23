# Plan: photometric calibration with Landolt standards (`calibrate` branch)

> Revised after three rounds of the user's answers, recorded at the end.
> Implementation waits for the user's go-ahead.

## Context

`pipnick_calibrate` will find the Landolt standards in a reduced data set and
measure the photometric transformation coefficients with
`pipnick/photometry/calibrate.py`. It needs the astrometry and photometry steps
to have run first. It will save the calibration in a form the user can edit,
assess it, and optionally apply it to the other observations of the night.

## What I checked

### The existing calibration code (`pipnick/photometry/calibrate.py`)

- **Model.** `m_inst - M = k0 + k1 X + k2 X C + sum_{n=1}^{order} k_{n+2} C^n`.
  The colour term and the airmass-colour term are always included. There is
  no way to drop a term, or to hold one fixed (e.g., a known extinction).
- **Fit.** A weighted linear least-squares solve (Harris et al. 1981). It
  returns no coefficient errors or covariance, and no residual statistics.
- **`photometric_calibration_coefficients_iterej`** needs an
  `(N_obs, N_band)` matrix, with every band measured for every
  "observation". That suits a matched multi-band star list. It does not suit
  frames taken one filter at a time, with a different number of standards in
  each. Fitting does not need it: the catalogue gives the colours, so each
  band can be fit on its own.
- **Bug, line 75.** `np.ones_like(mag_i.shape, dtype=bool)` makes a
  length-2 array, not a mask of `mag_i`'s shape. It is hidden only because
  callers always pass `gpm`.
- **`apply_photometric_calibration`.**
  - It starts the colour iteration from *instrumental* colours. This works,
    but the first guess is off by the zero-point difference.
  - The TODO at line 311 notes that it masks a whole row when any band is
    missing. Only the bands in each colour need to be valid.
- **Units.** The docstrings say `mag_inst` is in counts. The photometry
  tables give it per second. Nothing breaks, because the per-second
  convention is just absorbed into `k0`, but the docstrings should say so.

### The Landolt catalogues (queried from VizieR; copies in the scratchpad only)

| Catalogue | Stars | Coordinates | Proper motions | Size as ECSV |
|---|---|---|---|---|
| **II/183A**, Landolt (1992) | 526 | 1 s / 1″ precision | no | 73 kB |
| **J/AJ/137/4186**, Landolt (2009) | 595 | UCAC2, about 0.1″ | yes | 111 kB |
| **J/AJ/146/131**, Landolt (2013) | 349 | UCAC4 | yes (with errors) | 78 kB |

- The 2009 catalogue re-observes and supersedes the 1992 equatorial
  standards, with accurate coordinates and proper motions. The 2013
  catalogue adds fields at declination +50° and more equatorial stars. All
  three give V and the colours B−V, U−B, V−R, R−I and V−I, with errors and
  numbers of observations.
- **Matched against the real C2 2024-06-24 n1 photometry** (with a 6″
  radius):
  - Against 1992, the median offset is **3.7–4.3″** for the PG1530+057
    stars, because the coordinates are coarse and lack proper motions.
  - Against 2009, the offset is **0.03–0.35″**, and it recovers *more*
    stars: 3 per PG1530 frame instead of 2, and the 109_199 star in every
    frame where 1992 finds it in only one.
  - 2013 has no stars in these fields.
- **There are very few standards per frame.** There are 3 in each PG1530+057
  frame and 1 in each 109_199 and 109_231 frame, so about 5 distinct stars
  per band for the night.
  - The airmass covers only **1.185–1.356**.
  - The two 109_231 matches are 1.5–1.8″ off. The star is V ≈ 9.3, so it is
    probably saturated, and its centroid is biased.
  - Consequently, for this night, a model with both extinction and
    colour-airmass terms is barely constrained. This shapes the defaults
    below.

### Nominal extinction (from the Mt. Hamilton curve)

- **The curve.** `pypeit/data/extinction/mthamextinct.dat` has 36 rows,
  from 3200 to 10870 Å, in magnitudes per unit airmass. pypeit's
  `extinction_curves.txt` gives its site (lon −121.6428, lat 37.34139,
  elevation 1290 m) and warns "Generally grabbed from the Internet. Buyer
  beware". pypeit is BSD-3 licensed, so the file can be copied with
  attribution.
- **The bands.** I folded the curve through the speclite `bessell-UBVRI`
  responses (speclite 1.0.0, installed in the scratchpad only). The source
  is assumed to have a flat f_λ, and the response is photon-weighted:

  | Band | k (mag/airmass) | Curve covers the band |
  |---|---|---|
  | U | 0.604 | 84% (no data below 3200 Å) |
  | B | 0.300 | 100% |
  | V | 0.178 | 100% |
  | R | 0.120 | 100% |
  | I | 0.075 | 100% |

  These values are typical of Mt. Hamilton. U is the least certain.
- **Limits.** The effective extinction depends on the source spectrum. That
  dependence is what the `airmass_color` term models, and nominal values do
  not capture it. The Bessell curves also only approximate the Nickel
  filters.

## Proposed design

### 1. Packaged catalogue (`pipnick/data/catalogs/landolt.ecsv`)

- **Contents.** Merge Landolt (2009) and (2013) into one table:
  - `name`, `ra`, `dec` (deg, J2000 at epoch 2000.0), `pmra`, `pmdec`
    (mas/yr);
  - `U`, `B`, `V`, `R`, `I` and their errors, *derived* from V and the
    colours, with the errors propagated naively, i.e., ignoring the
    correlations between colours;
  - the original colours and their errors, the number of observations, and
    the source catalogue.
- **Duplicates.** A star in both catalogues keeps the 2013 values and notes
  this. The 1992 catalogue (II/183A) is not used.
- **Size and attribution.** About 200 kB. The table meta and the module
  docstring cite the papers and VizieR, and the CDS acknowledgement is
  included.
- **Building it.** A function, `pipnick.photometry.landolt.build_catalog()`,
  regenerates the file from VizieR. It is run by hand, not by the tests, and
  is marked `remote_data`.
- **Reading it.** `load_landolt()` reads the packaged file. Add
  `data/catalogs/*.ecsv` and `data/extinction/*.dat` to `package-data`.
- **Epoch.** Positions are propagated to each image's MJD with the existing
  `astrometry.services.propagate_proper_motion`. This matters little (the
  stars move about 1″ over 25 years at most), but it is cheap.
  - The catalogue therefore includes a `ref_epoch` column (2000.0).
  - When building it, check the units and convention of the Landolt proper
    motions: whether `pmRA` includes cos δ, and whether they are in mas/yr
    or arcsec/century. Convert them to what `propagate_proper_motion`
    expects.

### 1b. Packaged extinction curve and nominal coefficients

- **The file.** Copy `mthamextinct.dat` to
  `pipnick/data/extinction/mthamextinct.dat`, unchanged. Its attribution
  (pypeit, the site, and the "buyer beware" caveat) goes in the reading
  module's docstring and the docs.
- **Linking it to the camera.** Both camera TOMLs get
  `extinction = "mthamextinct.dat"` under `[observatory]`, and
  `NickelCamera` reads the curve from it.
- **`pipnick.photometry.extinction`:**
  - `load_extinction(file)` returns the curve;
  - `band_extinction(curve, band)` folds the curve through the speclite
    `bessell-<band>` response, photon-weighted, for a flat f_λ source. Where
    the band extends beyond the curve, the end values are held constant,
    and a warning gives the fraction of the band covered.
- **New dependency: `speclite`.** It is a small, pure-Python package (on
  PyPI; depends on numpy, astropy and pyyaml). It also has the SDSS
  responses, for later work with primed filters.
- **Using the values.** `--fix airmass=nominal` fixes the extinction terms at
  these values. An explicit `--fix airmass=B:0.28,...` still overrides them.
  When the user gets measured nominal coefficients from Lick, the camera TOML
  can list them per band, overriding the curve (`[observatory.extinction_k]`).
  This is an optional table, added only if the user wants it.

### 2. Finding the standards

1. **Load the photometry tables.** Use `reduction.load_setup` and
   `plan_products` for the science images, and their
   `photometry/<image>_phot.ecsv` tables.
   - A missing table is a warning: "run `pipnick_photometry`".
   - A table that is older than its image (`image_hash` differs from the
     image's `PNRAWHSH`/`PNTABH`) is an error, so that stale photometry
     cannot be used silently.
2. **Require celestial coordinates.** A table needs `wcs_method` set and
   `ra`/`dec` unmasked. Otherwise the image is skipped with a warning: "run
   `pipnick_astrometry`, then `pipnick_photometry`". If *no* table has
   coordinates, it is an error.
3. **Identify standards by matching, not by name.**
   - Every image with coordinates is matched against the catalogue, so a
     Landolt star in any field counts, whatever the `object` is.
   - Match radius: default 2″ (`--match_radius`). The pointing errors after
     the astrometry are at most about 0.3″.
4. **Select the measurements.** A matched measurement is used only if it
   passes all of the following:
   - `mag_inst` is unmasked;
   - `flags == 0`;
   - `flux_aper` is unmasked. It is masked if any pixel in the aperture is
     masked, and the reduction masks saturated pixels, so this rejects
     saturated stars. The rejection reason is "saturated or otherwise
     masked";
   - `group_size == 1`, i.e., isolated;
   - the catalogue error in the band is at most `--max_cat_err` (default
     0.05 mag);
   - the star has at least `--min_cat_obs` (default 2) Landolt observations.

   Each rejection is logged with its reason and kept in the residual table,
   masked.
5. **One calibration per setup and UT night.** The night is the UT date of
   each image's MJD; a Lick night falls within one UT date. A table that
   spans nights gives one calibration per night.
6. **`--images` limits only the standards used in the fit.** `--apply` acts
   on every science image of that setup and night.

### 3. Model and fit

- **Bands and colours.** The band map (Nickel filter to Landolt band) and
  the default colour for each band:

  | Filter | Band | Default colour |
  |---|---|---|
  | `U` | U | U−B |
  | `B` | B | B−V |
  | `V` | V | B−V |
  | `R` | R | V−R |
  | `I` | I | V−I |

  Each can be overridden (`--color B=B-V`). Any other filter (e.g., Andor
  `ip` and the other primed filters) is skipped with a warning: calibration
  currently requires UBVRI filters.
- **Which instrumental magnitude.** `--mag psf` (the default) uses
  `mag_inst`. `--mag aper` uses −2.5 log10(`flux_aper`/exptime) and logs a
  warning. The warning says that no aperture correction is applied, so the
  calibration is valid only for magnitudes measured in the same aperture
  (2 FWHM by default). The magnitude used is recorded in the calibration
  file, and `--apply` uses the same one.
- **Terms selected by name:** `zp`, `airmass`, `color`, `color2`,
  `airmass_color` (`--terms`).
  - **Default: `zp airmass color`.**
  - Any term can be *fixed* at a given value (`--fix airmass=B:0.25,V:0.15`,
    or from a file). The fixed terms are subtracted before the fit.
  - **Each term's source is recorded** as `"fit"`, `"fixed"` (by `--fix`) or
    `"fallback"`.
  - **Fixed and fallback terms carry no error.** The fitted errors, and the
    errors propagated from them, are statistical only. Fixing a term adds a
    systematic error that is not included. The docs, the TOML comments and
    the header comments say so.

  To support this, generalize `calibration_design_matrix` to take a list of
  terms.
- **Fit each band independently,** stacking every standard measurement in
  that band (star × image).
  - The weights are `1/(σ_inst² + σ_cat² + σ_floor²)`, where `σ_floor` is
    `--sys_err` (default 0.01 mag) and stops the brightest stars dominating.
  - Iterative sigma clipping (`--rej`, default 3σ, `--maxiter`), with each
    clipped point recorded.
  - Return the coefficients **and their covariance** (from the inverse of
    the normal matrix, scaled by the reduced χ² when it is above 1), the
    reduced χ², and the rms.
- **Sample checks**, per band, before and after rejection:

  | Check | Default | Outcome |
  |---|---|---|
  | Distinct stars < 3 (`--min_stars`) | 3 | Error, unless `--force` is set; then a warning. |
  | Measurements ≤ free terms | — | Always an error, even with `--force`, because the fit cannot be solved. |
  | Airmass span < `--min_airmass_span` while `airmass` or `airmass_color` is free | 0.2 | **Automatic fallback:** `airmass` is fixed at its nominal value, and `airmass_color` is dropped (fixed at 0); see below. `--no_fallback` makes this an error instead. |
  | Colour span < `--min_color_span` while a colour term is free | 0.3 mag | Error, suggesting fixing the term. A warning with `--force`. |
  | Fitted coefficient error > 0.1 mag/unit, or unphysical (e.g., negative extinction) | — | Warning. |

  - **An error affects only its band.** The band is skipped and recorded as
    failed, with the reason, in the TOML and the summary. The other bands
    go ahead.
  - **The run raises** only if no band is calibrated.

  On this n1 night, the airmass check would fire for every band (span
  0.17), so every band would fall back.
- **The extinction fallback must be impossible to miss.** A logged warning
  can easily be overlooked, so the fallback is recorded in every product the
  user reads:
  - **The calibration TOML:** the band's `source` entries record
    `"fallback"` for the affected terms. A `fallback_reason` gives why, e.g.,
    `"airmass span 0.17 < 0.2"`. A TOML comment above the band says that the
    extinction was *not measured*.
  - **The image headers:** `PNCALEXT` = `'fit'`, `'fixed'` or `'fallback'`.
  - **The residual table:** its meta records the same for each band.
  - **The QA figures:** a prominent "EXTINCTION NOT MEASURED — nominal
    value used" label on each affected band's panels.
  - **The terminal:** a warning when the fallback happens, and an
    end-of-run summary table (also in the log file) with an extinction
    column. Any fallback is repeated in a banner at the very end of the
    output.
  - **The docs:** a `.. warning::` admonition in `docs/calibrate.rst`
    explaining:
    - when the fallback happens and what it means;
    - that nominal values come from a site curve folded through generic
      Bessell filters, not from measurements of that night;
    - how to tell whether it happened (the keyword, the TOML field);
    - how to avoid it (standards over a wider airmass range) or turn it off
      (`--no_fallback`).

### 4. Assessment

- **`calibrate/<setup>_<night>_residuals.ecsv`.** One row per matched
  standard measurement:
  - image, star, band, airmass, catalogue magnitude and colour, `mag_inst`;
  - the calibrated magnitude, the residual and its error;
  - whether it was used, and the reason if it was not.
- **Per-band summary** in the log and the calibration file: the number of
  stars and measurements, the rms, the reduced χ², and the airmass and
  colour ranges.
- **QA figures** (`--qa`, matplotlib `Figure` as in `astrometry.qa`): per
  band, residuals against airmass, colour, catalogue magnitude and time,
  with the clipped points marked.
- **Leave-one-star-out residuals.** For each star, refit without any of its
  measurements, and predict them. The residual table gets a `resid_loo`
  column, and the summary gets the rms of the leave-one-out residuals.
  - Each leave-one-out refit applies the same fallback rule.
  - `resid_loo` is masked when the reduced fit has no remaining degrees of
    freedom (e.g., 3 stars and 2 free terms), or when it cannot be solved.
  - The summary gives the number of stars with a leave-one-out residual.
- **Drift (assessed, not fit).** Per band, and for all bands together,
  report the weighted linear slope of the residuals against time
  (mag/hour), with its error and significance. Warn when it is significant
  at more than 3σ and larger than `--drift_warn` (default 0.02 mag/hour). The
  QA figure shows the trend. A fitted time term is noted in the docs and
  code as an area for improvement.

### 5. Saved calibration (editable)

`rdxdir/calibrate/<setup>_<night>.toml`, written with a small writer for
this fixed structure, because the standard library writes no TOML:

```toml
[calibration]
setup = "..."
night = "2024-06-25"
catalog = "Landolt (2009, 2013)"
magnitude = "psf"          # which instrumental magnitude was calibrated
created = "..."
input_hash = "..."         # photometry table hashes + settings; unchanged inputs skip the refit
fit_hash = "..."           # hash of the fitted content below, to detect user edits

# WARNING: extinction NOT MEASURED for this band; nominal value used
[band.V]
filter = "V"
status = "ok"              # or "failed", with `reason`
color = "B-V"
terms = ["zp", "airmass", "color"]
coefficients = [ -24.50, 0.178, 0.03 ]
errors = [ 0.01, 0.0, 0.01 ]           # statistical only; fixed/fallback terms carry no error
source = ["fit", "fallback", "fit"]
fallback_reason = "airmass span 0.17 < 0.2"
covariance = [[...], [...], [...]]
n_stars = 5
n_meas = 12
rms = 0.021
rms_loo = 0.030
drift = [0.004, 0.006]     # mag/hour and its error; assessed, not applied
airmass_range = [1.185, 1.356]
color_range = [-0.2, 1.1]
images = ["PG1530+057_V_C014_2024-06-25", ...]
```

- **User edits.** A user can edit `coefficients` (or `terms`). To fix a
  term in a refit, use `--fix`.
- **Reruns:**
  - **Unchanged inputs:** if the `input_hash` is unchanged, the refit is
    skipped. The file, and the headers' `PNCALHSH`, are left alone.
  - **Changed inputs:** the night is refitted and the file rewritten.
  - **A hand-edited file** (its `fit_hash` no longer matches its contents)
    is kept, and is used by `--apply`. If its `input_hash` also no longer
    matches the inputs, a warning says the edits may hide changed
    photometry.
  - **`--overwrite`** always refits.

### 6. Applying the calibration (optional, `--apply`)

- **To the source tables.**
  - The calibrated magnitudes are written to *new* tables,
    `calibrate/<image>_cal.ecsv`, so that the photometry products are never
    altered.
  - The colour term needs each source's colour. Sources in the images of
    the same object and night are cross-matched by sky position (all bands
    at once). `apply_photometric_calibration` then iterates the colours,
    after its TODO is fixed and it starts from the colour corrected by the
    zero point.
  - **Repeated images in a band.** Each image is still calibrated and
    written separately. The colour used comes from the error-weighted mean
    of the source's calibrated magnitudes across images, in each band,
    iterated.
  - **`--default_color`** (per colour, default 0) is used for a source seen
    in only one band of its colour, or whose colour involves a band without
    a calibration. The source is flagged `color_assumed`.
  - Columns: `mag`, `mag_err` (statistical: from the instrumental error and
    the coefficient covariance), `color_used`, `color_assumed`.
- **Filters without a calibration.** Their images are skipped, for both the
  source tables and the headers, with a warning, and listed in the summary.
- **To the images: header keywords only.**
  - **No copy is needed.** The pixels are unchanged, so a full copy would
    duplicate the data for a few header cards. The astrometry already sets a
    precedent: it writes `PNAST*` into the reduced image headers in place
    (`fits.open(mode='update')`). The calibration does the same.
  - **Why this is safe:**
    - The photometry's rerun check (`_image_hash`) uses only `PNRAWHSH` and
      `PNTABH`, so the new keywords do not trigger re-measurement.
    - If the image is re-reduced, the keywords are lost. The TOML file is
      their record (as the sidecar is for the astrometry), so `--apply image`
      restores them.
    - Keywords are written only with `--apply image`. On such a rerun, they
      are rewritten only if `PNCALHSH` does not match the current
      calibration file.
  - **Keywords.** The model is `m_inst − M = k0 + k1·X + (k_c + k_xc·X)·C`.
    So the calibrated magnitude is
    `M = -2.5 log10(rate) + PNCALZP + PNCALCT·C`, where `rate` is in e-/s:
    - `PNCALZP = −(k0 + k1·X)`: the zero point at the image's airmass for
      C = 0;
    - `PNCALZE`: its statistical error (fixed terms contribute none);
    - `PNCALCT = −(k_c + k_xc·X)`: the colour coefficient at the image's
      airmass;
    - `PNCALCOL`: the colour, e.g., `'B-V'`;
    - `PNCALMAG`: `'psf'` or `'aper'`;
    - `PNCALSYS`: `'Vega'` (Landolt magnitudes are Vega-based);
    - `PNCALEXT`: `'fit'`, `'fixed'` or `'fallback'`: how the extinction was
      set;
    - `PNCALFIL`, `PNCALHSH`: the calibration file and its hash;
    - `PHOTZP`, a widely recognised name, equal to `PNCALZP`.

    The header comments say that `PHOTZP` assumes C = 0.

### 7. Script and pipeline

- **`pipnick/scripts/calibrate.py:Calibrate`**, entry point
  `pipnick_calibrate`. The options:
  - `-t` (required), `-r`;
  - `--images` (restrict the standards used);
  - `--terms`, `--fix` (including `airmass=nominal`), `--color`;
  - `--mag {psf,aper}` (default `psf`, matching `mag_inst`);
  - `--match_radius`, `--max_cat_err`, `--min_cat_obs`, `--sys_err`;
  - `--rej`, `--maxiter`;
  - `--min_stars`, `--min_airmass_span`, `--min_color_span`, `--force`,
    `--no_fallback`;
  - `--drift_warn`, `--qa`, `--overwrite`;
  - `--apply {catalog,image,both}`, `--default_color`.
- **`pipnick/pipelines/calibrate.py:calibrate_all`**, mirroring
  `photometry_all`.
- `pipnick/photometry/landolt.py` holds the catalogue I/O and matching.
- `pipnick/photometry/extinction.py` holds the extinction curve and the band
  coefficients.
- `pipnick/photometry/calibrate.py` is reworked freely:
  - named terms, fixed terms, and the covariance;
  - a per-band fit with rejection, replacing `_iterej`;
  - leave-one-out and drift helpers;
  - the `ones_like` bug, the per-colour masking TODO, and the docstrings
    fixed.

### 8. Tests (`pipnick/tests/test_calibrate.py`; no network)

- **Synthetic photometry tables** built from a handful of packaged Landolt
  stars, with known injected coefficients, noise and airmasses.
- **Recovery** of the coefficients within their errors, for each term set,
  with fixed terms honoured.
- **Checks:** each sample check raises or warns as specified, and
  `--force` downgrades the error to a warning.
- **Selection:** flagged, saturated, grouped and unmatched measurements are
  excluded with the right reasons.
- **Input errors:** missing coordinates and stale photometry tables give
  the documented errors and warnings.
- **The TOML file:** the round trip works; an edited file is detected and
  kept on rerun; `--overwrite` refits.
- **`--apply`:**
  - catalogue magnitudes are recovered for a synthetic two-band object;
  - the header keywords reproduce the catalogue magnitudes, and the pixels
    are unchanged;
  - a rerun with an unchanged calibration leaves the header alone;
  - the photometry still sees the image as up to date.
- **Leave-one-out and drift:** an injected outlier star stands out in
  `resid_loo`; an injected linear drift is detected and triggers the
  warning, and none is reported without one.
- **Extinction:** `band_extinction` reproduces the values above, and
  `airmass=nominal` fixes them. speclite's filter files ship with the
  package, so this needs no network.
- **Fallback:**
  - a narrow airmass span falls back to the nominal values, and drops a
    free `airmass_color`;
  - the fallback is recorded in the TOML, the header (`PNCALEXT`), the
    residual meta and the final summary;
  - fixed and fallback terms carry zero error;
  - `--no_fallback` raises an error instead.
- **Failures and reruns:**
  - a failing band is skipped while the others are calibrated, and the run
    raises only if every band fails;
  - unchanged inputs skip the refit and leave the file and headers alone;
  - a kept edited file with changed inputs gives a warning.
- **Leave-one-out with no degrees of freedom:** `resid_loo` is masked.
- **Aperture magnitudes:** `--mag aper` logs the warning that there is no
  aperture correction.
- **Other filters:** a non-UBVRI filter is skipped with its warning.
- **Regression tests:** the `ones_like` bug, and the per-colour masking
  (the TODO).
- **The packaged catalogue:** it loads, and its derived magnitudes are
  consistent with V and the colours.

### 9. Docs and release notes

- A new `docs/calibrate.rst`: the inputs, the catalogue and its citation,
  the model and terms, the checks, the outputs, editing, and applying. It
  also covers:
  - the extinction-fallback warning described in section 3;
  - the source and caveats of the Mt. Hamilton curve;
  - a note that the filter responses are speclite's generic Bessell UBVRI
    curves, to be replaced by Lick/Nickel-specific curves when available;
  - a note that zero-point drift is assessed but not fit, which is an area
    for improvement;
  - a note that the reported errors are statistical only: a fixed or
    fallback term adds a systematic error that is not included.
- API pages, and release notes in `docs/releases/0.1.0dev.rst`.

### 10. Real data (scratchpad only, FITS files only, no submissions)

`rdx_v2` is stale (the read noise changed), so first re-run the reduction,
astrometry and photometry on C2 n1 in the scratchpad. Astrometry may need
an astrometry.net submission, which needs your approval. Then calibrate
B, V, R and I with the defaults. Check that the extinction fallback is
reported everywhere specified in section 3, and look at the residuals, the
leave-one-out residuals and the drift.

## Decisions recorded (first round)

1. Ignore II/183A; package Landolt 2009 + 2013.
2. No aperture correction. `--mag aper` warns about this.
3. One calibration per setup and night.
4. Default terms: `zp airmass color`.
5. Nominal extinction: package `mthamextinct.dat` from pypeit, and fold it
   through the speclite filter responses. The user is also asking Lick for
   measured nominal coefficients.
6. Non-UBVRI filters are skipped, with a warning that calibration currently
   requires UBVRI.
7. Images get header keywords only. They are written in place, so no copy
   is needed.
8. Leave-one-out residuals are included.
9. Drift is assessed but not fit. It is noted as an area for improvement.
10. A hand-edited calibration file is kept on rerun unless `--overwrite` is
    set.
11. `calibrate.py` may be changed freely.

## Decisions recorded (second round)

A. **Automatic fallback** to the nominal extinction when the airmass span is
   too small. It must be made exceedingly clear: it is recorded in the
   products and the docs, not only in a logged warning (section 3).
   `--no_fallback` makes it an error.
B. **`speclite`** is a core dependency.
C. **The speclite Bessell filter curves** are used for now. The docs note
   that they can be replaced later by Lick/Nickel-specific curves.

## Decisions recorded (third round, after a full re-read)

- **A failing band** is skipped, and the run raises only if no band is
  calibrated.
- **`--force`** downgrades the minimum-stars and colour-span errors. It does
  not downgrade "measurements ≤ free terms".
- **One `source` per term** (`"fit"`, `"fixed"` or `"fallback"`) replaces
  `fixed` and `extinction_source`. Fixing a term for a refit is done only
  with `--fix`.
- **The fallback also drops a free `airmass_color`.**
- **No nominal error is added.** The errors are statistical; fixed terms add
  a systematic error, which is not modelled.
- **`input_hash`:** unchanged inputs skip the refit. A kept edited file with
  changed inputs gives a warning.
- **`--images`** limits only the standards. `--apply` covers every science
  image of the night whose filter has a calibration.
- **Repeated images:** each image is calibrated separately. The colours
  come from error-weighted means across images, and `--default_color` is
  used when a band of the colour is missing.
- **Leave-one-out:** `resid_loo` is masked when the reduced fit has no
  degrees of freedom.
- **Images in filters without a calibration** are skipped with a warning,
  and listed in the summary.


## Implementation progress (2026-10-23)

Implemented as planned: `pipnick_calibrate` (`pipnick/scripts/calibrate.py`),
`pipnick/pipelines/calibrate.py`, `pipnick/photometry/{landolt,extinction}.py`,
the rewritten `pipnick/photometry/calibrate.py`, the packaged catalog (942
stars, 206 kB) and extinction curve, `speclite` as a core dependency, 22 new
tests (`pipnick/tests/test_calibrate.py`), `docs/calibrate.rst`, the API
pages, release notes and `licenses/README.rst`.  Full suite: 145 passed,
1 skipped (also under `tox -e 3.13-test`); codestyle clean; no Sphinx
warnings.

Deviations from the plan:

- `--min_cat_obs` defaults to 1, not 2, and standards without a catalog error
  are used (with a warning giving their number).  The PG1530+057 standards
  have a single Landolt observation and no published errors, so the planned
  defaults would have rejected them all.
- `fit_calibration` raises if the fitted terms are degenerate (rank check),
  and the band fails with that reason.  Without it, a band observed at a
  single airmass gave an invalid fit (and NaN errors).
- The catalog values are rounded to 4 decimals (positions to 7) so that the
  single-precision VizieR values do not add spurious digits.

Real data (C2 2024-06-24 n1, re-reduced in the scratchpad as `rdx_v3`; no
astrometry.net submissions): runs end to end; the fallback is reported in
every product; reruns are no-ops.  Only R is calibrated (3 measurements,
rms 0.29 mag); B, V and I fail with fewer than 3 usable stars.  Causes:

- most standards have `group_size > 1` (rejected as blended), often from
  faint companions or spurious detections near bright stars;
- PG1530+057B has a bad column in its aperture;
- the 109_199 frames (those with stale header pointings) have fluxes that
  differ by factors of 2-4 between consecutive frames in the same filter and
  exposure time, which drives the color terms.

Open items for the user:

- whether to relax the `group_size == 1` criterion (e.g., reject only
  companions brighter than some fraction of the standard);
- the reduction and photometry rerun checks do not track camera
  configuration changes (e.g., the read noise), so `rdx_v3` needed
  `--overwrite` for both.

### Relaxing `group_size == 1` (2026-10-23)

Tested on C2 n1 (scripts and outputs in the scratchpad only).  Every source
grouped with a standard was detected in a later iteration of the PSF
photometry (2 or 3), within 0.2-1.4 FWHM, with 1-6% of the standard's flux:
residual detections from the PSF mismatch, not real neighbors.  Grouped
standards have aperture/PSF flux ratios of 1.02-1.05 (ungrouped: 1.00-1.01).
Repeat differences of PG1530+057 stars in back-to-back frames: PSF up to
0.045 mag; PSF with the companions merged, no better; aperture at most
0.012 mag.  The 109_199 frames differ by more than 1 mag in B and I with every
method, so they are a data problem.

Adopted (nominal): `group_size` is no longer used.  Residual detections (later
iteration, within `--residual_radius` = 1.5 FWHM) are counted (`n_residual`)
but do not exclude a standard; any other source within the aperture with more
than `--max_neighbor_flux` = 1% of the flux makes it `contaminated`.  On n1 this
keeps the same standards as ignoring the grouping.  All four bands now
calibrate, but with rms 0.2-0.9 mag and leave-one-out rms 0.5-1.5 mag, driven
by the 109_199 frames.  A dataset of better quality is needed for a robust
test.

Open: whether `--mag aper` should be the default; the PSF mismatch for bright
stars in the photometry.
