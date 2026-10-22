# Plan: camera configurations and a correct `pipnick_metadata` table

> Design record for the `metadata` branch. Status: approved; implementation
> proceeds phase by phase (see "Implementation phases"), and per-night
> verification results are recorded at the end of this file.

## Context

pipnick was developed on 2024 data from the **old Nickel camera** (Loral 2Kx2K,
UCAM controller; `obs/2024-06-24`, nights `n1`, `n2`). The 2026 commissioning
run used a **new Andor camera** (DZ936; `obs/2026-08-17`, nights `n1`–`n4`).
Detector properties are hard-coded in `pipnick/utils/nickel_data.py` and in two
camera classes in `pipnick/cameras/cameras.py`. The Andor class is a placeholder
copy of the old one and cannot even identify Andor files.

The datasets differ in almost every metadata convention:

| | Old camera (`NickelC2`) | New camera (`NickelAndor`) |
|---|---|---|
| Identify by | `INSTRUME='Nickel Direct Camera'`, `DSENSOR`, `CAMERAID=2` | `CAMID='headmodel DZ936…'`, `ANDORKIT`, `INSTRUME='Direct'` |
| Exposure time | `EXPTIME` (bias = 0) | `XPOSURE` (bias ≈ 0.2–0.5 s; no `EXPTIME`) |
| Filter | `FILTNAM` | `FILTER` |
| Binning | `CDELT1U/2U` = −2 | `CDELT1U/2U` = +2 (also 1, 4; `CBIN/RBIN` = log2) |
| Read speed | `READ-SPD` (80; 10 on failed reads) | `AMPCONF` ('1.0MHz' / '0.05MHz') |
| Raw shape (2x2) | 1056×1024; 32 overscan cols right; `DATASEC` present | 1044×1024; ~9 prescan + ~10 overscan cols, partial cols at both edges; no `DATASEC` |
| Frame type | inferable from `OBJECT` + `EXPTIME` | `OBSTYPE='DARK'` = shutter closed; otherwise free-text, stale `OBJECT` |
| Within a night | one setup | mixed binning (n3), mixed read speed (all), one window (n2039) |

The new camera's header content is **not final**; commissioning is ongoing. The
mixed setups within its nights reflect commissioning tests, and a typical night
will use a single setup. Even so, `pipnick_metadata` must identify and
categorize frames by setup. Section 7 recommends header changes for the new
camera that would simplify supporting both cameras. The design reads those
keywords whenever they are present, so header improvements need only config
edits, not code changes.

**Goals for this branch:**

1. Move all camera properties into per-camera **TOML** configuration files,
   selected automatically from the raw headers, or explicitly by camera name.
   Include support for an archival bad-pixel mask.
2. Make `pipnick_metadata` build a correct, typed metadata table for each of
   the six nights. Each night is processed independently.

**Scope:** configuration plus metadata. Reduction is rewired only to *read*
geometry, the BPM, gain/read noise and saturation from the configuration.

**Deferred:**
- Measuring the Andor gain, read noise and bad pixels; placeholders are flagged
  in the config.
- Running full end-to-end reductions.

**Decisions made with the user:**
- Configs are TOML (stdlib `tomllib`, no new dependency).
- Frame classification is chosen on the command line: header-only, or header
  plus image statistics.
- Users exclude frames by **commenting out rows** with `#`. This is verified to
  work with astropy's ECSV reader, and column types are preserved.

## Design

### 1. Camera configuration files: `pipnick/data/cameras/{nickel_c2,nickel_andor}.toml`

They are packaged automatically: `include-package-data` plus setuptools_scm
include git-tracked files. They are loaded with `importlib.resources`. All
geometry is in **unbinned** pixels; binned sections are derived. Sections:

- `[camera]`:
  - `name`, `description`;
  - `pixel_size` (µm) and `plate_scale` (arcsec per unbinned pixel: C2 0.185;
    Andor ≈0.1665, from the header CD matrices);
  - `detector_shape`, `saturation`, `fwhm_guess`.
- `[identify]`: header keyword → regex. All must match.
  - C2: `INSTRUME = '^Nickel Direct Camera$'`, `DSENSOR = '^Loral'`.
  - Andor: `CAMID = '^headmodel DZ936'`.
- `[keywords]`: maps generic names to an **ordered list** of header keywords;
  the first one present wins. Generic names: object, frametype, obstype,
  exptime, filter, airmass, ra, dec, date, xbin, ybin, readspeed, gain,
  readnoise, datasec, biassec, window, temperature.
  - Example for Andor: `exptime = ["EXPTIME", "XPOSURE"]`,
    `datasec = ["DATASEC"]`.
  - If the new camera's headers adopt the standard keywords in section 7, they
    take precedence automatically. The geometry, readout and frame-type rules
    below remain the fallback when those keywords are absent.
- `[geometry]`: unbinned prescan, data and overscan column ranges (and rows),
  plus partial-pixel handling at boundaries. C2 also falls back to the header
  `DATASEC`/`COVER`.
- `[[readout]]` table, one entry per read-speed value. Fields: `mode`
  (`slow`/`fast`), `gain`, `readnoise`, `placeholder = true/false`.
  - C2: `READ-SPD=80` → existing 1.8 e⁻/ADU and 10.7 e⁻.
  - Andor: `1.0MHz` → fast, `0.05MHz` → slow. Values are placeholders; a
    warning is logged when they are used.
- `[bpm]`:
  - `columns`, `regions` (polygons) and optional `file` (archival unbinned BPM
    image), all in unbinned coordinates.
  - C2: the existing bad columns and corner polygons, converted from 2x2
    to unbinned. Add column 1 (always saturated), the partly illuminated row 1,
    and the edge roll-off, after confirming them in the data.
  - Andor: partial and pre/overscan boundary columns only for now. The `file`
    hook is in place for a measured mask later.
- `[frametypes]`: an ordered list of rules (`[[frametypes.rule]]`). The first
  match wins. Each rule can test: `obstype` regex, exposure-time bounds, and a
  case-insensitive `object` regex.
  - If the header provides an explicit frame type (the `frametype` keyword
    list, e.g. `IMAGETYP`), it is mapped through `[frametypes.map]` and used
    before any rule.
  - C2:
    - `exptime == 0` → bias;
    - `OBJECT` `^flat$`/`sky flat` → skyflat;
    - `dome flat` → domeflat;
    - `^focus` → focus;
    - `^point` → pointing;
    - `test` → test;
    - otherwise → science.
  - Andor:
    - `OBSTYPE=DARK` and `exptime < 1` → bias;
    - `OBSTYPE=DARK` → dark;
    - `dome ?flat` → domeflat;
    - `twilight|sky ?flat` → skyflat;
    - `focus` → focus;
    - `pointing` → pointing;
    - `test` → test;
    - otherwise → science.
- `[quality]`: thresholds used to auto-comment frames:
  - maximum detector temperature;
  - bias-level range;
  - saturated or empty fractions (used only in stats mode).

### 2. Camera classes: rewrite `pipnick/cameras/cameras.py`

- One `NickelCamera` class constructed from a config. Remove the
  `NickelC2`/`NickelAndor` subclasses and the `eval(f'cameras.{name}')` calls.
- A registry, `available_cameras()` / `get_camera(name)`, built by discovering
  the TOML files.
- `identify_camera(hdr)` returns a `NickelCamera`. It raises if zero or more
  than one camera matches.
- Methods:
  - `parse_metadata(hdr)` returns typed, normalized values;
  - `binning(hdr)` returns positive ints;
  - `readout(hdr)` returns mode, gain and read noise;
  - `sections(hdr)` returns binned data and overscan FITS sections, including
    windows via `CRVAL1U/2U` and `WINDOW`;
  - `bpm(hdr)` gives the unbinned mask, OR-binned to the frame's binning and
    cut to its window;
  - `classify(meta, stats=None)`.
- Reuse `pipnick/utils/masking.py` (`create_mask`, `mask_polygons`) for the BPM.
- Keep `get_header`.

### 3. Metadata table: move `pipnick/utils/dir_nav.py` to **`pipnick/metadata.py`**

The module is really the metadata-table layer, so the name `dir_nav`
("directory navigation") no longer fits. The new name follows PypeIt's
`pypeit/metadata.py`, and the module sits at the package level because the
reduction and photometry pipelines both depend on it.

- **Moved:** `build_camera_datatable`, `build_metadata`, `frame_path`,
  `default_table_name` and `default_reduced_file_name`, plus the new
  `read_metadata`/`write_metadata`.
- **Moved elsewhere:** `unzip_directories`, which builds `Fits_Simple` objects
  for `photometry/fit.py` and `utils/display_fits.py`, goes to
  `pipnick/utils/fits_class.py`.
- **Deleted:** `comment_out_rows`, `norm_str` and `create_exclusion_func`, which
  are unused (checked with grep). The new comment-preserving I/O replaces
  `comment_out_rows`.
- **Callers to update:** `pipelines/reduction.py`, `pipelines/photometry.py`,
  `scripts/metadata.py`. Alternatives considered: `pipnick/utils/metadata.py`
  and `pipnick/rawdata.py`.

The functions:

- `build_camera_datatable`:
  - identifies the camera **per file**;
  - raises a clear error on mixed cameras, unless `--camera` is given, in which
    case non-matching files are skipped with a warning;
  - fails cleanly if no files are found;
  - builds **typed** columns. The current `np.asarray` turns every column into a
    string and `None` into `'None'`.
- Columns (the same for both cameras):
  - from the header: `file`, `frametype`, `object`, `filter`, `exptime`,
    `airmass`, `ra`, `dec`, `date` (DATE-BEG) and `mjd`;
  - setup: `xbin`, `ybin`, `readmode`, `naxis1`, `naxis2`, `window`;
  - `setup`: a short label (`A`, `B`, …, as in PypeIt), assigned to each
    distinct combination of camera, binning, read mode and window. Frames are
    identified and categorized by setup even though a typical night has only
    one. Each label's definition is stored in `meta['setups']`, e.g.
    `A: {xbin: 2, ybin: 2, readmode: slow, window: full}`. A bias or dark can
    only calibrate frames in the same setup; filter is handled separately when
    flats are matched;
  - `note`: why a row was auto-commented or reclassified.
- `meta`: `camera`, `rawdir`, `classify` (mode used) and the pipnick version.
- Classification modes:
  - `--classify header` (default) applies the config rules and header quality
    checks, e.g. the warm Andor frames in n1.
  - `--classify stats` also reads each image and computes the median and the
    saturated and zero fractions over the data section. Then:
    - empty frames (all zeros: C2 d1061, n2 d1030/1031) and saturated test
      frames are auto-commented;
    - a closed-shutter frame at bias level is reclassified as bias, whatever its
      `OBJECT` says (e.g. Andor n1032–1041, labelled 'test');
    - open-shutter frames whose level contradicts their type get a `note`.
- **Comment-preserving I/O**, replacing the legacy `comment_out_rows`:
  - `write_metadata(tbl, path, commented=mask)` writes the ECSV, then prefixes
    `# ` to the rows to be commented out.
  - `read_metadata(path, include_commented=False)` reads it back. With
    `include_commented=True` it un-comments the rows (identified by file name)
    to recover them.
  - Re-running with `--append` keeps existing rows *as written*: user edits to
    `frametype` and similar columns, and commented-out rows. Only new files are
    added. `--overwrite` rebuilds from scratch.
- Remove unused imports (`pandas`, `focus_label`).

### 4. `pipnick_metadata` CLI: `pipnick/scripts/metadata.py`

- Arguments:
  - `-d/--rawdir` (a single night's raw directory);
  - `-t/--raw_table`;
  - `--camera NAME` (overrides auto-detection; choices come from the registry);
  - `--ext` (default `.fits`);
  - `--classify {header,stats}`;
  - `-o/--overwrite`, `-a/--append`.
- Fix the copied module docstring.
- Log a summary so users can sanity-check a night at a glance:
  - a frames-per-setup × frametype table, and each setup's definition;
  - the number of commented-out rows;
  - warnings for setups that have science frames but no matching bias or flat.

### 5. Retire `pipnick/utils/nickel_data.py` and the unused `pipnick/utils/nickel_masks.npz`

Point their consumers at the camera config:

- `pipelines/astrometry.py`: plate scale, FOV and FWHM guess, RA/DEC keys, bad
  columns. Its separate import error, `organize_files`, stays out of scope and
  is noted.
- `photometry/fit.py`: plate scale.
- `utils/graphs.py`: use the image shape instead of `ccd_shape`.
- `utils/fits_class.py`: filter and exptime keywords come from the camera. Drop
  the broken `get_masks_from_file` fallback.
- `photometry/psf_photometry.py`: remove the unused `bad_columns` import.
- Delete `docs/api/pipnick.utils.nickel_data.rst`; `make apirst` regenerates the
  API pages.

### 6. Minimal reduction rewiring: `pipnick/pipelines/reduction.py`

- `init_ccddata`/`process_frame` take overscan and trim sections, the BPM,
  saturation, gain and read noise from the camera object for each frame's
  binning and read mode. This replaces `oscansec`/`COVER`/`DATASEC` and fixed
  shapes.
- Science frames are `frametype == 'science'`; the `'None'` sentinel is gone.
- Guard: if the uncommented rows contain more than one `setup` (relevant for
  Andor n3 and the n2039 window), raise an error naming the setups. The user
  then comments out frames, or reduction grouping is extended in a later
  branch.

### 7. Recommended header changes for the new camera

The goal is self-describing headers that use the same standard keywords as the
old camera, or standard FITS/IRAF conventions. pipnick then needs fewer
camera-specific rules. The old camera's headers are archival and won't change.
In priority order:

1. **Explicit frame type.** Add `IMAGETYP` (or make `OBSTYPE` reliable), set by
   the observing software from the observer's selection: `bias`, `dark`,
   `domeflat`, `skyflat`, `object`, `focus`, `pointing`, `test`.
   - This is the biggest single improvement. Today, frame type must be inferred
     from shutter state plus free-text `OBJECT`, which is often stale, e.g.
     biases labelled 'test' and flats labelled from a previous exposure.
2. **`EXPTIME`** in seconds, as the requested exposure time, with `0` for
   biases. Keep `XPOSURE` as the measured value.
3. **Detector sections:** `DATASEC`, `BIASSEC` and `TRIMSEC`, plus `PRESEC`
   for the prescan. The header's own COMMENT cards say these are missing. Also
   `DETSEC`/`CCDSEC` (or `LTV1/LTV2`, `LTM1_1/LTM2_2`) so binned and windowed
   frames map back to detector coordinates.
4. **Binning** as `CCDSUM = '2 2'` (an IRAF convention) or positive
   `XBINNING`/`YBINNING`. The `CDELT1U/2U` sign differs between cameras, and
   `CBIN`/`RBIN` are log2 values.
5. **`GAIN` and `RDNOISE`** in e⁻/ADU and e⁻ for the readout mode actually
   used. Note that the old camera's `GAIN` is a DCS index, not e⁻/ADU.
6. **Readout mode as a named keyword**, e.g. `READMODE = 'slow'/'fast'`,
   alongside `AMPCONF`.
7. **Camera/detector identification:**
   - `DETECTOR` (e.g. 'Andor DZ936 SN31382') and a stable `INSTRUME`;
   - keep `TELESCOP = 'Nickel'`;
   - adding `TELESCOP` and `DETECTOR` to the old camera's data is not possible,
     but new data should have them.
8. **`FILTNAM` vs `FILTER`:** pick one name. If the observatory keeps `FILTER`,
   writing `FILTNAM` too would let legacy tools keep working. Also record the
   filter-wheel position (cf. `FILTRAW` on the old camera).
9. **Times:** add `DATE-OBS` and `MJD-OBS` (standard) alongside `DATE-BEG` and
   `DATE-END`.
10. **`SHUTTER`** state (open/closed/auto), and **lamp state that reflects
    reality**: `LAMPS` is always 'Off', even for dome flats.
11. **Software hygiene:**
    - clear or reset `OBJECT` when the exposure type changes;
    - reset the `OBSNUM` counter at the start of each night (the n2170 file in
      n3 was taken before the reset);
    - don't write a file when a readout fails.

The plan **does not depend on any of these changes**. Each new keyword is added
to the front of the corresponding `[keywords]` list, so it's used when present,
and the current inference rules remain the fallback for existing data.

### 8. Tests (focused on the new code) and docs

Add `pipnick/tests/` cases using **synthetic headers** (no data files):
- config loading and registry;
- identification of each camera, and the unknown and ambiguous cases;
- keyword normalization (exptime, filter, binning sign, read mode);
- frame-type rules for representative `OBJECT`/`OBSTYPE`/exptime combinations
  taken from both datasets;
- binned and windowed section math against the measured Andor and C2
  geometries;
- BPM binning;
- setup assignment, with one setup and with several;
- keyword priority: a synthetic Andor header with the recommended standard
  keywords (`EXPTIME`, `DATASEC`, `IMAGETYP`, `CCDSUM`) uses them over the
  fallbacks;
- an ECSV comment and append round trip.

Docs:
- add release notes to `docs/releases/0.1.0dev.rst`;
- document the TOML schema in the `pipnick/cameras` module docstring.

## Critical files

- **New:**
  - `pipnick/data/cameras/nickel_c2.toml`, `nickel_andor.toml`
  - `claude/metadata_plan.md`
  - `pipnick/tests/test_cameras.py`, `test_metadata.py`
- **Rewritten:** `pipnick/cameras/cameras.py`, `pipnick/scripts/metadata.py`
- **Moved and rewritten:** `pipnick/utils/dir_nav.py` → `pipnick/metadata.py`
  (`unzip_directories` → `pipnick/utils/fits_class.py`)
- **Updated:** `pipnick/pipelines/reduction.py`, `pipnick/pipelines/astrometry.py`,
  `pipnick/photometry/fit.py`, `pipnick/utils/fits_class.py`,
  `pipnick/utils/graphs.py`, `pipnick/tests/test_mask.py` (camera-config
  values)
- **Deleted:** `pipnick/utils/nickel_data.py`, `pipnick/utils/nickel_masks.npz`
- **Reused:** `pipnick/utils/masking.py` (`create_mask`, `mask_polygons`),
  `frame_path` and `default_table_name` (moving to `pipnick/metadata.py`),
  `scriptbase.ScriptBase.init_log`

## Implementation phases

0. Save this plan as `claude/metadata_plan.md`. (Done.)
1. TOML configs, plus the camera class, registry and identification. Add tests.
2. Move `dir_nav.py` to `pipnick/metadata.py` and update its callers. Then
   add typed metadata extraction and normalization, frame-type rules, and setup
   identification. Add tests.
3. Comment-preserving ECSV I/O, `--append`/`--overwrite`, `--classify stats`,
   and the CLI. Add tests.
4. Verify on all six nights (below), iterate on the rules and thresholds, and
   record the per-night results in `claude/metadata_plan.md`.
5. Retire `nickel_data.py` and the npz, and update their consumers.
6. Minimal reduction rewiring, the setup guard, docs, release notes and
   codestyle.

## Verification

- `tox -e 3.13-test,codestyle` passes.
- Run `pipnick_metadata` on each night, writing tables to a scratch directory,
  never into the data directories:
  - `obs/2024-06-24/n{1,2}/raw`
  - `obs/2026-08-17/n{1,2,3,4}`
- **Nothing from the data directories is added to this repository.** Only
  the FITS files are read; the observing logs and other notes are not, because
  they contain credentials.
- Check each night in both `--classify header` and `--classify stats` mode:
  - camera auto-detected; column types correct; nothing is `'None'`;
  - frame-type and setup counts are reported per night, and the user checks
    them against the observing logs (or provides redacted copies);
  - known cases handled:
    - C2 n2 biases with `OBSTYPE='OBJECT'`;
    - C2 n2 d1012 ('Bias' at 10 s → not bias);
    - C2 empty frames (d1061, n2 d1030/1031) commented out in stats mode;
    - Andor n1032–1041 (bias labelled 'test') → bias in stats mode;
    - Andor n2132–2135 (OBJECT filter ≠ FILTER) use the `FILTER` value;
    - every night gets setups whose definitions match the logs; the C2 nights
      and the typical Andor nights have one or two setups;
    - Andor n3 has a separate setup for each binning × read speed, and the
      stray n2170 lands in the 1x1 slow setup;
    - n2039 window gets its own setup;
    - Andor n1 warm frames are flagged.
- Append test: comment out a row, change one `frametype`, then rerun with
  `--append`. The edits survive and no rows are duplicated.
- Smoke test: `pipnick_reduce -t <C2 n1 table>` still runs through bias and flat
  creation with geometry and BPM from the config, and the setup guard triggers
  on Andor n3.
