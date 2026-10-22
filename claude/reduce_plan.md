# Plan: rework `pipnick_reduce` around combination groups (`reduce` branch)

> Design record for the `reduce` branch.  Status: implemented and verified;
> awaiting review.  Progress is recorded at the end of this file.

## Context

`pipnick_metadata --split` now writes one metadata file per setup, with a
`comb_id` column marking which frames to combine. The reduction pipeline in
`pipnick/pipelines/reduction.py` (`reduce_all`) predates this:

- It groups frames itself: all biases; flats by filter, preferring sky over
  dome.
- It ignores darks.
- It processes science frames one at a time, writing `{stem}_rdx.fits`.
- It decides whether to redo an output only by whether the file exists.

**Goals of this branch:**

1. Use `comb_id` to decide what is combined.
2. Enforce rules that make calibrations deterministic:
   - one superbias;
   - one dark per exposure time;
   - one flat per filter for each flat type;
   - a sky flat is preferred over a dome flat for the same filter.
3. Name outputs predictably, e.g. `M11_C002_2026-08-19.fits` for science.
4. Record in each output's header the raw frames it was made from.
5. Reprocess an existing output only when its inputs have changed, so
   reductions can be rerun as data come in.

**Decisions made with the user:**

- **Science groups:** combined without registration. The frames must share a
  filter. Different exposure times are allowed, with a warning; frames are
  scaled to the median exposure time before combining.
- **Darks:** master darks are made, one per exposure time, but not yet applied
  to science frames.
- **Science file names** use the **UT date** of the group's first frame.
- **One setup per reduction directory.** Calibration names stay as they are
  (`Bias.fits`, `SkyFlat_B.fits`). Writing over, or reusing, a product from
  another setup is an error.
- **Reprocessing trigger:** an output is redone if its raw frames changed, or
  if any calibration it used was rebuilt from different raw frames.
- **Input:** `pipnick_reduce` requires a per-setup file with `comb_id`. For a
  single-file table, it stops with an error pointing to
  `pipnick_metadata --split`.
- **`--append` is removed from `pipnick_reduce`.** New frames are added with
  `pipnick_metadata --split --append`, and `pipnick_reduce` is simply rerun.

## Design (`pipnick/pipelines/reduction.py` unless noted)

### 1. Input and grouping rules

- **Signature:** `reduce_all(table, rdxdir=None, overwrite=False)`. The table is
  required; `rawdir` and `append` are removed. It reads the table with
  `build_metadata(filename=table)` and then `check_single_setup`, both
  existing functions.
- **New `check_comb_groups(metadata, metafile)`.** It considers uncommented
  rows only. It collects every violation and raises one
  `PipnickDataModelError` listing them all:
  - the `comb_id` column is missing (point to `--split`);
  - a group mixes frame types;
  - there is more than one bias group;
  - a dark group mixes exposure times (rounded to 0.1 s), or two dark groups
    share an exposure time;
  - a dome or sky flat group mixes filters, or two groups of the same flat type
    share a filter;
  - a science group mixes filters.
- **Warning, not an error:** a science group whose exposure times differ.
- **Frame types that are not processed:** focus, pointing, test and unknown
  frames are ignored if they're uncommented, with an info message. Only bias,
  dark, flat and science groups are processed.

### 2. Output names

- **Calibrations:** `Bias.fits`, `Dark_<exptime>s.fits` (e.g. `Dark_120s.fits`,
  `Dark_4.3s.fits`; exposure time rounded to 0.1 s and formatted with
  `{:g}`), `SkyFlat_<filter>.fits`, `DomeFlat_<filter>.fits`.
- **Science:** `<object>_C<comb_id:03d>_<UT date>.fits`.
  - `<object>` is the most common `OBJECT` in the group; ties go to the
    earliest frame.
  - `<UT date>` is the date part of `date` for the group's earliest frame (by
    `mjd`).
- **Name sanitizing** (helper `_safe_name`), for the object and the filter:
  spaces become `_`, and any character other than letters, digits, `.`, `+`
  and `-` is dropped. For example, `NGC 7331 I` becomes `NGC_7331_I` and `i'`
  becomes `i`.

### 3. Provenance in headers

Every product, calibration or science, gets:

- **HISTORY lines:** one per raw frame, `pipnick raw frame: n3101.fits`. Science
  outputs also get one per calibration used, e.g.
  `pipnick calibration: Bias.fits`.
- **Keywords for the reprocessing check:**
  - `PNSETUP`: the per-setup tag, e.g. `2x2_fast_full`;
  - `PNCOMBID`: the group number;
  - `PNNRAW`: the number of raw frames;
  - `PNRAWHSH`: a hash of the sorted raw file names (SHA-1, first 16 hex
    digits).
- **Calibration keywords on science outputs** (and on flats and darks, which
  depend on the bias): `PNBIAS`/`PNBIASH` and `PNFLAT`/`PNFLATH`, giving each
  calibration's file name and `PNRAWHSH`. The value is `none` when no
  calibration was available.

### 4. Reprocessing logic

New `_needs_update(ofile, setup_tag, raw_hash, calib_hashes, overwrite)`:

- Returns `False` only if the file exists, `overwrite` is not set, and the
  file's `PNRAWHSH` and calibration hashes all match.
- If the file's `PNSETUP` differs from the current setup, it raises
  `PipnickDataModelError` (one setup per directory). Files without `PNSETUP`,
  i.e. older outputs, are treated as needing an update.
- The log says whether each output was "up to date" or "reprocessed", and why.

Order of processing:

1. Bias first.
2. Darks and flats next. Each depends on the bias's hash, so a rebuilt bias
   cascades to them.
3. Science last. Each depends on the bias hash and the hash of the flat chosen
   for its filter: the sky flat if one exists, otherwise the dome flat,
   otherwise none, with a warning.

**Stale outputs:** after processing, any `*.fits` in the reduction directory
with this setup's `PNSETUP` that wasn't produced or confirmed in this run (for
example, a group that has since been merged into another) is reported in a
warning. Stale files are never deleted.

### 5. Processing and combination

- **Keep** `init_ccddata` and `process_frame` as they are.
- **Refactor `stack_frames`** into `process_frames`, which returns the processed
  frames, and `combine_frames(frames, scale=None)`, which does the existing
  sigma-clipped mean with optional per-frame scale factors.
  - Biases are combined unscaled.
  - Flats keep the existing scaling by the mean.
- **Darks:** bias-subtracted and combined; not applied to science frames.
- **Science groups:** each frame is processed (bias, flat, cosmic rays).
  - A single frame is written as it is.
  - Multiple frames are scaled to the median exposure time and combined;
    `EXPTIME` is set to the median.
  - **Pointing check:** warn if any frame's header RA/Dec is more than 3
    binned pixels from the first frame's, using the camera plate scale times
    the binning, since frames are not registered.
- **Return value:** the list of science output paths, whether written or
  already up to date. The docstring already says this; the code currently
  returns the metadata table.

### 6. Command line and other code

- **`pipnick/scripts/reduce.py`:** `-t/--raw_table` is required; `-r/--rdxdir`
  and `-o/--overwrite` stay; `-d/--rawdir` and `-a/--append` are removed. Help
  text describes the per-setup input and the rerun behavior.
- **Removed** from `reduction.py`: `get_super_bias`, `get_super_flats` and
  `stack_frames`, replaced by group-based functions:
  - `_reduce_bias`, `_reduce_darks`, `_reduce_flats`, `_reduce_science`;
  - helpers for names, hashes and header provenance.
- `metadata.default_reduced_file_name` is no longer used by the reduction. It
  stays for `pipelines/photometry.py`, which needs updating for the new names
  in a later step. That pipeline is already mid-refactor and broken.
- **Docs:** update the `reduce_all` docstring and add release notes in
  `docs/releases/0.1.0dev.rst` (breaking: `pipnick_reduce` options and output
  names).
- **Code style:** follow `AGENTS.md`, without black.

### 7. Tests (`pipnick/tests/test_reduction.py`)

Use synthetic 4x4-binned Andor frames (512×522) written with
`write_andor_image` from `test_metadata.py`, and per-setup files built with
`build_metadata(split=True)`. Every assertion has a message.

- **`check_comb_groups`:** each violation raises; a mixed-exposure science
  group only warns.
- **Name helpers:** `_safe_name`, the dark name formatting, and the science
  name (most common object, UT date of the first frame, `C%03d`).
- **End to end:** biases, V dome flats, V sky flats, an R dome flat, a 120 s
  dark pair, and a science group of two V frames:
  - the expected products are written;
  - the science output uses the sky flat;
  - HISTORY lists the raw frames, and the keywords are set.
- **Rerun:** produces no reprocessing; file modification times are unchanged.
- **Adding a frame to a science group** (by editing `comb_id`) reprocesses only
  that output.
- **Adding a bias** reprocesses the bias, the flats, the dark and the science
  output.
- **A product from another setup** in the directory raises an error.
- **Stale outputs** are reported but not deleted.
- **Single-file table:** fails with the `--split` message.

## Critical files

- `pipnick/pipelines/reduction.py`: the main rework.
- `pipnick/scripts/reduce.py`: command-line options.
- `pipnick/tests/test_reduction.py`
- `docs/releases/0.1.0dev.rst`
- **Reused:**
  - from `pipnick/metadata.py`: `build_metadata`, `read_metadata`,
    `frame_path`, `update_metadata_meta`, `setup_file_name` (or the setup tag
    logic behind it);
  - from `pipnick/pipelines/reduction.py`: `check_single_setup`,
    `init_ccddata`, `process_frame`;
  - from `pipnick/cameras/cameras.py`: the camera `plate_scale` and
    `binning`, for the pointing check.

## Verification

- `tox -e 3.13-test,codestyle` passes.
- **Real data,** scratch outputs only, FITS-only reading. Reduce the old camera
  n1 per-setup file (`n1_2x2_slow_full.ecsv`) into a scratch directory:
  - **Products:** `Bias.fits`; `SkyFlat_{B,V,R,I}.fits`; 26 science files named
    `<object>_C<id>_2024-06-25.fits`, each with HISTORY and keywords.
  - **Rerun:** nothing is reprocessed; file modification times are unchanged.
  - **Combining:** edit `comb_id` to merge two same-filter frames of one target
    and rerun. One combined output is written; the now-unused single-frame
    output is reported as stale; the others are untouched.
  - **Wrong setup:** reducing another setup's file into the same directory
    fails with the one-setup error.

## Progress

- **Implemented** as planned in `pipnick/pipelines/reduction.py`,
  `pipnick/scripts/reduce.py`, `pipnick/tests/test_reduction.py` and
  `docs/releases/0.1.0dev.rst`. Small differences from the plan:
  - `_reduce_bias`, `_reduce_darks` and `_reduce_flats` are one function,
    `_reduce_calibration`; science groups use `_reduce_science`.
  - Darks and flats also get `pipnick calibration:` HISTORY lines for the bias.
  - `metadata.setup_tag` was added (and is used by `setup_file_name`) to build
    the `PNSETUP` tag.
  - The one-setup check scans every FITS file in the reduction directory before
    any processing, so a conflict stops the run before anything is written.
  - Two groups that would write the same file name are an error.
  - Follow-up: in flat names, a primed filter is written with a `p` (`i'`
    becomes `ip`), so it stays distinct from `i`; `_filter_name` does this.
    Object names still just drop quotes.
- **Tests:** `tox -e 3.13-test,codestyle` passes (77 tests).
- **Real data** (2024-06-24 n1, `n1_2x2_slow_full.ecsv`, scratch outputs only):
  - Products: `Bias.fits`, `SkyFlat_{B,V,R,I}.fits`, and 26 science files
    `<object>_C<id>_2024-06-25.fits`; about 30 s.
  - Rerun: all 31 up to date, no modification times changed.
  - Merging d1033+d1034, and d1057–d1060: only the two outputs were
    reprocessed. The mixed exposure time and pointing offset warnings were
    issued, and four stale outputs were reported.
  - Reducing an Andor setup into the same directory raised the one-setup error
    and wrote nothing.
- Follow-up: science names include the filter,
  `<object>_<filter>_C<comb_id:03d>_<UT date>.fits` (e.g.
  `M11_V_C002_2026-08-19.fits`), using `_filter_name`.
- Follow-up: the table is the authority for the object, filter, exposure time
  and coordinates (`_table_values`, `_apply_table_values`). These are written
  to the output headers, with differing raw values kept in `ORIGOBJ`,
  `ORIGFILT`, `ORIGEXP`, `ORIGRA` and `ORIGDEC`. Their hash (`PNTABH`) is part
  of the reprocessing check.
