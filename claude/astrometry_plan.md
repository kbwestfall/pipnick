# Assessment and plan: update the astrometry pipeline (`astrometry` branch)

> Design record for the `astrometry` branch. Status: implemented and tested
> with synthetic data; awaiting a real astrometry.net run (needs the user's
> permission). Progress is recorded at the end of this file.

## Context

`pipnick_astrometry` solves the WCS of the reduced images with the
astrometry.net web service (nova.astrometry.net) through
`astroquery.astrometry_net.AstrometryNet`. Each solve is a remote job on
their servers. The code predates the metadata and reduce branches and no
longer runs.

**Constraints from the user:**
- **Ask the user for permission before submitting any job to astrometry.net
  during development.** This covers every run of the pipeline that could
  submit a job, however small (one image, `--max_submissions 1`, a retry,
  `--resolve`, or `--upload_image`). Before each such run, state:
  - which images could be submitted, and how many;
  - the command;
  - the expected number of jobs, from a `--dry_run` beforehand.

  Wait for explicit approval. An approval covers only the run it was given
  for, never later runs. Without approval, use `--dry_run`, the cache, or
  bootstrapping only.
- Submit to astrometry.net as rarely as possible during development.
- No test may submit a job.

## Assessment of the current code

### `pipnick/pipelines/astrometry.py`

- **Does not import.** It imports `organize_files` from `pipnick.metadata`,
  which no longer exists. `pipnick_astrometry --help` still works (the import
  is inside `main`), but every real run fails.
- **Expects the old directory layout:** `maindir/reduced/<obj_dir>/*_red.fits`.
  It writes zeroed copies of each image to `processing/mod-astro-input/` and
  solutions to `astrometric/wcs/` (`astrometric/corr/` is created but never
  used).

  The reduce branch instead writes one flat `rdxdir` per setup, named
  `<object>_<filter>_C###_<date>.fits`. Each file has `MASK` and `UNCERT`
  extensions and `PN*` provenance keywords. Its object, filter and RA/Dec
  come from the metadata table.
- **Frame selection uses the retired exclusion options**
  (`use_table`, `excl_files`, `excl_objs`, `excl_filts`). In the script these
  are `type=list`, so a string becomes a list of characters.
- **Uploads data as well as positions.** `solve_from_image` is used:
  - In astroquery 0.4.11 it still detects sources locally with photutils and
    submits a source list, but that path is deprecated. Once it is removed,
    the whole image is uploaded.
  - Separately, the pipeline writes masked copies of every image to disk
    first (`make_cleaned_inputs`).
- **Bugs:**
  - In the retry loop over timed-out submissions, `output_path` is whatever
    the last image set. Every late solution is saved to the same, wrong file.
  - `get_astrometric_solves` uses `output_path` before assigning it.
- **Partly reusable:**
  - The search hints:
    - the scale bounds come from `camera.plate_scale` × binning;
    - the field width comes from `NAXIS1`;
    - the center comes from the header RA/Dec. Sexagesimal is fine, since
      astroquery parses hour/deg, and these values now come from the table.
  - `publicly_visible='n'`.
  - The two-pass timeout/`monitor_submission` idea.
- **Rerun logic:** a solution is reused only if a file with that name exists.
  There is no link to the image's contents.

### `pipnick/scripts/astrometry_script.py`

- The API key is a positional argument, so it lands in the shell history and
  in `ps` output.
- The other options are the retired ones above.

### Related code

- **`pipnick/pipelines/final_calib.py`** (marked incomplete) reads
  `maindir/astrometric/`.
  - It keys solutions by `name.split('_')[0]`, which breaks with the new names
    and with object names that contain underscores.
  - It has no script entry point.
- **`pipnick/pipelines/photometry.py`** is still on the old naming (known,
  deferred).
- **`pipnick/photometry/starfind.py::detect_sources`** (photutils
  `DAOStarFinder`, isolation and Kron metrics) is available for local source
  detection. It may be heavier than astrometry needs.
- **Shared code to reuse:**
  - `reduction.check_single_setup`, `reduction._plan_products`,
    `reduction._group_indices`, `metadata.build_metadata`,
    `metadata.setup_tag` and `cameras.get_camera`;
  - the `PN*` provenance pattern (`_needs_update` style).
- **Network guard already in place:** `pytest-remotedata` comes with
  `pytest-astropy` (in the `test` extra).
  - Without `--remote-data`, which `tox.ini` never passes, it blocks
    internet access from tests.
  - It also skips tests marked `remote_data`.
- **No tests exist** for astrometry.

## Proposed steps

1. **Inputs follow the reduction.**
   - New signature: `astrometry_all(table, rdxdir=None, resolve=False, ...)`,
     reading the same per-setup metadata file as `reduce_all`.
   - The science outputs are found with the same naming logic, by making
     `reduction._plan_products` (or a science-only wrapper) public.
   - Images that don't exist yet are skipped with a warning.
   - The `PNSETUP` of each image is checked.

2. **Isolate the remote service behind a small solver interface.**
   - A class wraps `AstrometryNet`, with `solve(sources, width, height,
     hints)` returning a header or a submission ID on timeout, and
     `monitor(id)`.
   - Everything else is local and testable: detection, hints, header updates
     and bookkeeping.
   - Tests pass in a fake solver.

3. **Detect sources locally and submit only a source list**
   (`solve_from_source_list`), never the image.
   - Use the reduced image's `MASK` extension instead of writing zeroed
     copies.
   - Keep the N brightest isolated sources.
   - This sends less data, gives faster and more reproducible solves, and
     avoids the deprecated astroquery path.

4. **Hints.** The center comes from the table RA/Dec in the header. The radius
   and scale bounds come from the camera plate scale, binning and image size,
   as now.

5. **Validate each solution locally** before accepting it:
   - the pixel scale is within bounds;
   - the solved center is within the hint radius of the header pointing;
   - most detected sources map inside the field.

   Rejected solutions are logged and not applied.

6. **Store solutions and decide when to resubmit.** This is the key to
   minimizing submissions. Proposed:
   - Cache each returned solution as a small sidecar header file. It records
     the hash of the submitted source list and the hints, plus the
     astrometry.net submission and job IDs.
   - Apply the WCS to the reduced image header, with a `PNASTRO*` provenance
     keyword (status, job ID, cache hash).
   - On rerun:
     - an image whose WCS is already applied and current is skipped;
     - an image rewritten by `pipnick_reduce` (header replaced) has the
       cached solution re-applied, *without* resubmitting, if its source list
       is unchanged;
     - a new submission happens only when there is no cached solution for the
       current source list, or with `--resolve`.

7. **Robust retry handling.** Submit all images first. Then monitor the
   timed-out submissions, each mapped to its own image (fixing the
   `output_path` bug). Record failures, so that they aren't retried every run
   unless asked.

8. **API key.** Read it from an environment variable (e.g.
   `ASTROMETRY_NET_API_KEY`) or astroquery's own config, with an optional
   `--api_key` flag. Never take it as a positional argument, and never log it.

9. **Command line** (`pipnick/scripts/astrometry.py`, renamed from
   `astrometry_script.py`):
   - `-t` (per-setup table, required), `-r` (rdxdir) and `--resolve`;
   - `--upload_image`: upload a masked copy of each image instead of a source
     list;
   - `--dry_run`: detect sources and compute hints, write what would be
     submitted, and submit nothing. This is for checking inputs on real data
     during development.

10. **Testing without submissions.**
    - Synthetic reduced images with stars placed through a known TAN WCS, in
      the reduced-image format (`MASK` extension, `PN*` keywords).
    - A fake solver returns a header built from the known WCS, or raises a
      timeout and later returns from `monitor`, or fails.
    - The tests cover:
      - detection;
      - the hints;
      - validation (a good solution, and a wrong scale or center);
      - applying the WCS and provenance keywords;
      - cache reuse after the image is rewritten, with no new "submission";
      - `--resolve`;
      - the retry mapping;
      - the dry run;
      - `--upload_image`: the fake solver receives a masked copy, and the
        temporary file is removed;
      - that the API key is never logged.
    - A test asserts that the default solver is never constructed during
      tests, as a backstop alongside `pytest-remotedata`.
    - No `remote_data` tests are added.

11. **Development on real data.** Use `--dry_run` first. **Before any run
    that could submit a job, ask the user for permission** (see the
    constraints above), with `--max_submissions` set to the approved number.
    Rely on the cache and bootstrapping from then on. Outputs go to the scratchpad only, reading
    FITS files only.

12. **Docs.** Add release notes, covering the new options and the script
    rename. `final_calib.py` is deferred.

## Addendum A: solutions from matched points

**Goal:** a direct way to make a solution. A user (or the pipeline) supplies
pixel coordinates and celestial coordinates for a set of sources, and gets a
WCS back.

**What exists:** `astropy.wcs.utils.fit_wcs_from_points(xy, world_coords, ...)`
already fits a TAN (or other projection) WCS to matched points, with optional
SIP distortion, using `scipy.optimize.least_squares`. It is suitable; nothing
needs porting from astrometry.net, whose `fit-wcs` and `tweak` do the same job
in C.

**A version trap:**
- On astropy main (#20383, merged 2026-10-16), `fit_wcs_from_points` has an
  `origin` argument. Omitting it emits a `FutureWarning` and assumes 1-based
  pixels.
- The installed astropy 8.0.2 has no `origin`, and the fit is internally
  0-based. Testing it here showed:
  - 0-based input, evaluated with origin 0, gives an exact round trip;
  - 1-based input is off by 1.4 pixels.

**Plan:**
- **`fit_wcs(x, y, coords, sip_degree=None, sigma=None)`** in
  `pipelines/astrometry.py`:
  - It takes 0-based pixel coordinates (numpy convention, as photutils
    returns them).
  - It passes `origin=0` when the installed astropy accepts it (checked once
    with `inspect.signature`), and calls it without otherwise. That is
    correct on both versions.
  - It rejects outliers iteratively and refits.
  - It returns the WCS and the fit statistics: number of points used, RMS
    residual in arcsec, and the rejected points.
  - It requires a minimum number of points: 3 for TAN, more with SIP.
- **`apply_wcs(image, wcs, method, ...)`** writes a solution to a reduced
  image with the same provenance keywords as any other solution.
  - `PNASTMTH` gives the method: `astrometry.net`, `points`, or
    `bootstrap`.
  - `PNASTRMS` and `PNASTNPT` give the fit quality.
  - A solution from points is never replaced by a rerun unless `--resolve` is
    set.
- **Python only for now** (user decision). There is no command-line option.
  `fit_wcs` and `apply_wcs` are the interface, and `fit_wcs` is also the core
  of the bootstrapping below.
- **Tests:** points generated from a known WCS (with and without noise and
  outliers) are recovered. The test should pass on both the installed astropy
  and a version with `origin`.

## Addendum B: bootstrapping solutions between frames

**Premise:**
- Nickel's camera orientation and plate scale are fixed for a setup, and the
  required accuracy is a few tenths of an arcsec.
- Frames of the same field over a night differ mainly by a shift (pointing),
  plus at most a small rotation.
- So one astrometry.net solution can serve many frames.

**Proposed strategies, in order of preference:**

1. **Group by pointing first.**
   - Cluster the science images of a setup by their header pointing, using
     the table RA/Dec. Images whose pointings are within a radius (`--group_radius`,
     default 1 arcmin) form one pointing group.
   - In each group, at most one image is submitted: the reference, preferring
     the image with the most detected sources.
   - This alone reduces submissions to about one per field per night.

2. **Frame-to-frame transfer within a group** (local only, no catalog):
   - Detect sources in the new image (the same detection used for
     submissions).
   - Match them to the reference image's sources in pixel space:
     - Find the shift by voting on all pairwise offsets: take the peak of a
       2D histogram of `dx`, `dy` within the maximum expected offset, set by
       the pointing difference plus a margin.
     - Match by nearest neighbour within a tolerance.
     - Fit a rigid transform (shift plus rotation, scale fixed at 1) with
       outlier rejection.
     - Rotation is expected to be about 0, so the fit is well constrained
       with as few as the default minimum of 4 matches.
   - Map the matched reference pixels to the sky through the reference
     solution, then call `fit_wcs` (Addendum A) on the new image's pixels.
     With few matches (fewer than about 6), fit only a shift and rotation
     with the reference scale held fixed, instead of a full CD matrix.
     This carries over the reference's scale, rotation and any distortion
     without extra modelling.

3. **Transfer across groups, using the pointing error.**
   - For a group with no solved image, build a predicted solution: the CD
     matrix (scale and rotation) of a solved image in the same setup, with the
     center taken from the header pointing.
   - Correct it with the telescope's pointing offset, i.e. the difference
     between the solved and header centers. This is measured from the
     night's solved images (median, recorded in the log). It is usually
     stable over a night, and better near the solved images.
   - If the new field overlaps a solved field, refine with strategy 2.
     Otherwise, see strategy 4.

4. **Catalog refinement with Gaia** (confirmed by the user):
   - Query Gaia through `astroquery.gaia` for the field: a cone around the
     predicted center, a little larger than the image, limited to stars
     brighter than a magnitude cut and to the brightest few hundred.
   - Propagate the positions to the observation epoch with the Gaia proper
     motions. This is cheap, and keeps fast-moving stars from spoiling
     matches at the 0.1 arcsec level.
   - Starting from the predicted solution of strategy 3, match detections to
     Gaia sources by shift voting and nearest neighbour, then `fit_wcs`.
   - Cache each query as an ECSV file in the reduction directory, keyed by
     the cone center (rounded) and radius. Nearby images reuse a cached cone
     that covers them, so each field is queried once per night.
   - A Gaia query is a remote request but not an astrometry.net job. The
     query sits behind a catalog interface like the solver. Tests use a fake
     catalog built from the synthetic truth, and never query Gaia
     (`pytest-remotedata` blocks it as well).
   - A field that can't be matched to Gaia falls back to astrometry.net.

**Acceptance and fallback:**
- A bootstrapped solution is accepted only if:
  - at least 4 sources are matched (default, configurable with
    `--min_matches`; may be raised later);
  - the RMS residual is below 0.3 arcsec (default, configurable with
    `--max_rms`);
  - the scale and rotation agree with the reference to within small
    tolerances.
- Otherwise the image falls back to the next strategy, and finally to an
  astrometry.net submission.

**Provenance and reruns:**
- A bootstrapped image records `PNASTMTH = 'bootstrap'`, the reference
  image, and the reference solution's hash.
- If the reference is re-solved (or its cached solution changes), the
  dependent images are redone on the next run, like the calibration cascade
  in `pipnick_reduce`.
- Bootstrapping uses only local data, so redoing it costs nothing.

**Order in `astrometry_all`:**
1. Apply cached and points solutions.
2. Form the pointing groups.
3. For each group without a solution, choose a reference:
   - a solved image within the radius if there is one;
   - otherwise submit one image.
4. Bootstrap the rest of the group.
5. Retry the remaining failures through the fallback chain.

`--max_submissions` caps the number of jobs per run. It is useful during
development, e.g. `--max_submissions 1`. `--dry_run` reports the grouping and
which images *would* be submitted.

**Tests**, all with the fake solver and fake catalog, and no network:
- Synthetic images of one field, offset by up to about 1 arcmin and rotated
  slightly, form one pointing group. Exactly one fake submission is made, and
  each bootstrapped WCS matches the truth to < 0.1 pixel.
- The shift-voting matcher recovers known offsets, including with spurious
  sources.
- A failed bootstrap (too few matches) falls back to a submission.
- A field with no overlap is solved from the fake Gaia catalog with no
  submission, and the cached cone is reused for a nearby image without a
  second query.
- Proper-motion propagation moves a star by the expected amount.
- Re-solving a reference cascades to its dependents.
- `--max_submissions` and the dry-run report are respected.

## Decisions (confirmed by the user)

- **Where the solution lives:** applied to the reduced image header, with a
  sidecar cache of each solution (step 6).
- **Submission mode:** a source list is the default (step 3). The user can
  override with `--upload_image`:
  - A temporary copy of the image, with masked pixels set to 0, is uploaded
    with `force_image_upload=True`, and astrometry.net detects the sources.
  - The copy is deleted afterwards.
  - The cache key for these submissions is a hash of the uploaded data and
    the hints.
- **Script rename:** `pipnick/scripts/astrometry_script.py` becomes
  `pipnick/scripts/astrometry.py`, to match the other script names. Update the
  `pipnick_astrometry` entry point in `pyproject.toml`. `photometry_script.py`
  is renamed later, not on this branch.
- **`final_calib.py` is deferred.** Step 12 is only the release notes.
- **Solutions from points:** Python functions only (`fit_wcs`, `apply_wcs`)
  for now; no command-line option.
- **Bootstrapping defaults:** a 1 arcmin group radius, at least 4 matched
  sources (may be raised later) and an RMS residual of at most 0.3 arcsec.
- **Gaia:** query Gaia through astroquery for catalog refinement, with
  cached results, and never in tests.

## Critical files

- `pipnick/pipelines/astrometry.py` (rewrite);
- `pipnick/scripts/astrometry.py`: renamed from `astrometry_script.py` (with
  `git mv`), plus the entry point in `pyproject.toml`;
- `pipnick/pipelines/reduction.py` (expose the science-product listing);
- `pipnick/tests/test_astrometry.py` (new; synthetic data in
  `pipnick/data/tests/` only if a file is needed);
- `docs/releases/0.1.0dev.rst`.

## Verification

- `tox -e 3.13-test,codestyle` passes. No test reaches the network:
  `pytest-remotedata` blocks it, and the fake-solver backstop catches the rest.
- On real data, `--dry_run` on the C2 2024-06-24 n1 reduction lists the
  sources and hints for each image, with no submissions.
- With approval, one real run on a small subset (one or two images):
  - solutions are applied and validated;
  - a rerun makes zero submissions;
  - rerunning after `pipnick_reduce --overwrite` re-applies cached solutions
    with zero submissions.

## Progress

- **Implemented:**
  - `pipnick/astrometry/` (new package): `detect.py`, `fit.py` (`fit_wcs`,
    the astropy-version shim, `predict_wcs`, `wcs_from_transform`), `match.py`
    (shift voting, nearest-neighbour matching, rigid fit) and `services.py`
    (`AstrometryNetSolver`, `GaiaCatalog`, proper-motion propagation).
  - `pipnick/pipelines/astrometry.py` rewritten (`astrometry_all`,
    `apply_wcs`, `group_by_pointing`).
  - `reduction.load_setup` and `reduction.plan_products` are shared with the
    reduction.
  - Script renamed to `pipnick/scripts/astrometry.py`; entry point and API
    docs updated.
- **Differences from the plan:**
  - The algorithms live in a new `pipnick/astrometry` package (mirroring
    `pipnick/photometry`) instead of all in the pipeline module.
  - Gaia-solved images are marked `PNASTMTH = 'gaia'`.
  - While a submission is pending, no other image is submitted, since the
    pending solution may make further submissions unnecessary.
  - Pending submissions are cached and monitored on the next run, rather than
    resubmitted.
  - `astroquery.gaia` contacts the Gaia archive on import, so it is imported
    only when a query is made.
- **Tests:** `pipnick/tests/conftest.py` blocks every astrometry.net and Gaia
  method in all tests. `pipnick/tests/test_astrometry.py` has 13 tests using
  synthetic star fields, a fake solver and a fake catalog. `tox -e
  3.13-test,codestyle` passes (94 tests).
- **Real data, dry run only:** C2 2024-06-24 n1 (scratchpad reduction) has 22
  science images in 4 pointing groups. At least 1 submission is needed. No
  job has been submitted.

## Coordinate frames (design note)

**Problem:** the first real run crashed after the submission.
`SkyCoord.spherical_offsets_to` refused to compare the telescope pointing
(ICRS) with the center of the astrometry.net solution (FK5).

**Findings:**
- astrometry.net writes only `EQUINOX = 2000.0`, with no `RADESYS`
  (`util/sip_qfits.c:134` in the astrometry.net source).
- Its API has no setting for the output frame (see the astroquery
  `show_allowed_settings` list).
- By the FITS WCS standard, astropy reads such a header as FK5 J2000, and
  writes `RADESYS = 'FK5'` when it saves the header again.
- The solutions are fit to index files built from catalogs on the ICRS
  (Tycho-2, 2MASS, Gaia), so they are effectively ICRS. FK5 J2000 and the
  ICRS differ by at most tens of mas.

**Decisions (user-approved):**
1. Relabel astrometry.net solutions as ICRS: `services.icrs_header` sets
   `RADESYS = 'ICRS'` and removes `EQUINOX` and `RADECSYS`. The pipeline
   applies it to new results (`_handle_result`) and to cached results
   (`_restore_solution`); the cache written before this fix is labeled FK5.
2. `fit.image_center` always returns ICRS coordinates (a real conversion, not
   a relabel). This guards against user-provided solutions in other frames.

**Documentation:**
- module docstrings in `services.py` and `pipelines/astrometry.py`;
- the `icrs_header`, `image_center`, `AstrometryNetSolver` and `apply_wcs`
  docstrings;
- `docs/astrometry.rst` (rewritten, with a "Coordinate frames" section);
- the release notes.

**Tests:**
- The fake solver now returns headers as astrometry.net does (`EQUINOX`
  only).
- `test_icrs_header` and `test_image_center_frame` check the relabelling and
  the conversion.

## Real astrometry.net submissions (each approved by the user)

1. 2024-06-24 C2 n1, `--max_submissions 1`: submitted the 100-source list for
   `109_231_I_C020_2024-06-25.fits`. Solved (SIP order 2) and cached. The run
   then crashed in bootstrapping because of the frame mismatch above (fixed).
2. Rerun with `--max_submissions 0` (no approval needed), after the frame
   fixes. A second frame problem was found and fixed first: `SkyCoord.icrs`
   keeps the FK5 `equinox`, which astropy treats as a different frame, so
   `fit.to_icrs` builds a clean ICRS coordinate. An astrometry.net solution
   already applied with the FK5 label is relabeled from the cache.
   - 14 of 22 images solved with 0 submissions: 1 from the astrometry.net
     cache, 10 bootstrap, 3 Gaia. RMS residuals 0.02-0.08 arcsec.
   - The 8 unsolved frames (109_199, C024-C031; raw d1045-d1052) have a stale
     header pointing (the same as the 109_231 frames). An offline wide-search
     match to Gaia (Gaia query only) puts them at 17:45:02.6 -00:29:22, 5.4
     arcmin from the header pointing, with 0.05-0.09 arcsec fits. The
     pipeline correctly rejected the poor matches it found within its 30
     arcsec search.
   - Options: the user corrects the RA/Dec in the metadata table (supported,
     and propagated by `pipnick_reduce`), or the pipeline adds a wide-search
     fallback with stricter acceptance. To be decided.
3. Option 2 (user-approved): a wide-search Gaia fallback for incorrect
   pointings.
   - `_bootstrap_gaia(wide=True)`: search of `wide_search` arcmin (default 10),
     at least `wide_min_matches` sources (default 10), `PNASTPOE` and a
     warning.
   - Images with the same stale pointing (within 1 arcsec) reuse the measured
     error (`_predicted_center`). Wide-solved images are excluded from the
     pointing-offset median.
   - Two fixes found on real data:
     - Bootstrapping now uses the best overlapping solved image
       (`_nearest_solved`, preferring images that were not themselves
       bootstrapped), not the fixed group reference. Narrow-overlap fits
       (4-10 matches) had given 0.2-0.36 arcsec full-field errors.
     - `_invalidate_chains` propagates the invalidation along chains of
       bootstrapped solutions.
   - Result (`--max_submissions 0`): 22 of 22 solved (1 astrometry.net, 17
     bootstrap, 4 Gaia). The stale frames used 7-71 matches, RMS 0.07-0.12
     arcsec.
4. Option 1 verified on a copy: corrected the RA/Dec of d1045-d1052 in the
   table. `pipnick_reduce` reprocessed only those 8 images, and
   `pipnick_astrometry --wide_search 0 --max_submissions 0` solved all 22
   with the usual matching. Documented in `docs/astrometry.rst` ("Incorrect
   telescope pointings").

## Imports and products of the services (user request)

**Imports:**
- `import astroquery.gaia` creates a module-level `GaiaClass()`, which makes
  one HTTP request to the Gaia archive for status messages. That costs about
  0.7 s per import, prints a banner, and cannot be disabled.
- `GaiaCatalog` now uses `astroquery.utils.tap.core.TapPlus` (the client
  `GaiaClass` is built on) with the same server, imported at the top of
  `services.py`. Constructing it makes no request.

**Network blocking in tests:** `pytest-remotedata` blocks the internet only
with `remote_data_strict = true`, which was missing. It is now in
`pyproject.toml`. The earlier claim that it blocked internet access in tests
was wrong; `conftest.py` was the only guard. `conftest.py` now also blocks
`TapPlus.launch_job`.

**Products are kept as returned:**
- Gaia results are saved as returned (`gaia_<n>.vot.gz`; TapPlus keeps the
  `.gz` name), with the index `gaia_queries.ecsv`.
- astrometry.net solutions are saved as returned (`anet_<hash>.fits`, via
  `Header.tofile`), with the index `astrometry_net.ecsv` holding the status,
  submission ID, solution file and date.
- The relabelling, proper-motion propagation and column selection are applied
  on use.
- Tests compare the saved files with what was returned, byte for byte or card
  for card.

**Fixes found on real data while verifying:**
- **Normal search radius:** NGC 3982 B has a 31 arcsec pointing error and
  only about 8 Gaia matches. The 30 arcsec search was too tight, so the
  default `search` is now 60 arcsec (`--search`).
- **Pointing errors from any method:** images bootstrapped from a
  stale-pointing image were not flagged, and their 5.4 arcmin offsets skewed
  the median pointing offset. `_pointing_error` now flags any solution more
  than `search` from its table pointing: it writes `PNASTPOE`, warns, and
  excludes the image from the offset model. All 8 stale frames are now
  reported, and there are no false warnings.
- **Catalog margin:** the margin around the predicted image is applied only in
  the wide search.

**Result:** both real-data copies solve 22 of 22 with 0 submissions.

## Iteration tools and logging (user request)

- **Troubleshooting section** in `docs/astrometry.rst`. The "one good
  solution" bullet was clarified: one solved image per setup usually
  suffices (other groups are matched to Gaia), and solving one image per
  failing pointing group helps when that fails.
- **Options:**
  - `--images` (file names or comb_ids; `--resolve` applies only to them);
  - `--qa` (figures in `astrometry/qa/`, made with matplotlib `Figure`, not
    pyplot);
  - `--max_sources`, `--detect_thresh`, `--fwhm`.
  - `--gaia_mag_limit` was skipped (user decision).
  - The summary is written on every run (no `--summary` flag).
- **Logging:**
  - numeric failure messages;
  - two recorded reasons per image (`reason` for local solving,
    `submit_reason` for submission), combined in the summary with
    suggestions;
  - the end-of-run summary table, also in `astrometry/summary.ecsv`;
  - the median pointing-offset report;
  - a warning for fewer than 10 detections;
  - debug-level match details (`match_frames` now returns `nvotes`).
- **Real data** (`--max_submissions 0`): 22 of 22 solved and 21 QA figures.
  The figures show spurious detections (galaxy knots in NGC 3982, edge and
  masked-column artifacts), which supports the deferred source vetting.
