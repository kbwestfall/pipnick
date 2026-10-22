# Plan: one metadata file per instrument setup (`pipnick_metadata --split`)

> Design record for the `setups` branch.  Status: implemented; see "Progress"
> at the end of this file.

## Context

`pipnick_metadata` writes one ECSV table for a night. Each frame is labeled
with a setup (`A`, `B`, …), defined by camera, binning, readout mode and readout
window, and the definitions are stored in `meta['setups']`
(`pipnick/metadata.py`: `assign_setups`, `build_metadata`).

`pipnick_reduce` can only reduce one setup at a time: its setup guard in
`pipnick/pipelines/reduction.py` stops on a multi-setup table. A night with
several setups therefore has to be split by hand, by commenting out rows.

The goal is a second output mode that writes **one file per setup**. Every row
in such a file shares the same read speed, binning and window, so those values,
plus the image size they imply, move out of the columns and into the file's
ECSV metadata. The current single-file mode stays the default.

**Decisions made with the user:**

- **File names** describe the setup:
  `<base>_<xbin>x<ybin>_<readmode>_<window>.ecsv`. The window part is `full`
  for the full readout, and `win00`, `win01`, … for distinct custom windows,
  numbered in order of first appearance. Example: `n2_2x2_fast_full.ecsv`,
  `n2_2x2_fast_win00.ecsv`.
- **Moved to metadata:** the setup label and definition (`xbin`, `ybin`,
  `readmode`, `window`, plus the window name), and `naxis1`/`naxis2`.
- **Reading back** restores those values as columns, so `pipnick_reduce` and
  the rest of the code see the same table layout from either mode.
- **`--append`** is supported in the per-setup mode.

## Design (all in `pipnick/metadata.py` unless noted)

### File layout

A per-setup file has:
- the usual columns, minus `setup`, `xbin`, `ybin`, `readmode`, `window`,
  `naxis1` and `naxis2`;
- the usual `meta` (`camera`, `rawdir`, `classify`, `rdxdir`), but not the
  multi-setup `setups` dictionary;
- a `meta['setup']` dictionary instead. Example:
  `{label: B, camera: NickelAndor, xbin: 2, ybin: 2, readmode: fast,
  window: '[960:1205,840:1069]', window_name: win00, naxis1: 123, naxis2: 115}`.

Commented-out rows work exactly as in a single-file table. A setup whose frames
are all rejected (e.g., the old camera's failed readouts) still gets a file,
with every row commented out. Otherwise `--append` would treat those frames as
new each time and add them again.

### New and changed functions

- `setup_file_keys = setup_keys + ['naxis1', 'naxis2']`: the columns that move
  to `meta['setup']`. These are the quantities that define a setup, plus the
  image size, which they determine.
- `assign_setups(tbl, camera_name, setups=None)` (existing): also assigns each
  setup a `window_name`. It is `'full'` for the full readout. Each distinct
  custom window gets the next unused `win<NN>` (two-digit, zero-padded,
  starting at `win00`), and names already in `setups` are reused. Setups that
  differ only in binning or read mode but share the same window share its
  name.
  Labels and window names therefore stay stable on append.
- `setup_file_name(base, definition)` → `Path`. Builds
  `<stem>_<xbin>x<ybin>_<readmode>_<window_name><suffix>` from the base name's
  stem and suffix. If no base is given, the default is
  `default_table_name(camera, 'raw')`.
- `split_by_setup(tbl, commented)`:
  - yields `(definition, subtable, subcommented)` for each setup;
  - drops `setup` and the `setup_file_keys` columns;
  - sets `meta['setup']`, including `naxis1`/`naxis2`;
  - raises if a setup's image size is not uniform, which would indicate a bug.
- `read_metadata(filename, include_commented=False, expand=True)`: if the file
  has `meta['setup']`, it restores the dropped columns at their usual positions
  (order from `metadata_columns`, types preserved). It also rebuilds
  `meta['setups']` as `{label: definition}`, so the reduction setup guard and
  `summarize` work unchanged. `expand=False` returns the compact file layout.
- `update_metadata_meta` reads with `expand=False`, so updating a per-setup
  file (e.g., the `rdxdir` that `pipnick_reduce` adds) keeps the compact
  layout.
- `find_setup_files(base)`: returns the existing per-setup files for a base
  name. It uses `<stem>_*<suffix>`, keeping only files whose `meta` has a
  `setup` entry, so unrelated files with similar names are ignored.
- `build_metadata(..., split=False)`: refactored into three steps.
  1. **Obtain the combined table and its commented-row flags.**
     - *Build:* as now.
     - *Append, single-file:* as now.
     - *Append, split:* read every file from `find_setup_files` with
       `include_commented=True, expand=True`, stack them into one table, and
       rebuild `meta['setups']` from each file's `meta['setup']`. Then follow
       the existing append path: add only new files, keep edits and commented
       rows, keep existing labels and window names.
  2. **Write:** either one file (as now) or, with `split=True`, one file per
     setup from `split_by_setup`, named by `setup_file_name`.
     - In append mode, every per-setup file is rewritten from the combined
       table, with values and commented rows preserved, and new setups get new
       files.
     - With `overwrite`, existing per-setup files that aren't rewritten (stale
       setups) are listed in a warning, not deleted.
  3. **Summarize:** `summarize` on the combined table, logged once as now. In
     split mode it also lists each file written and its frame counts.
  - **Return value in split mode:** `(list of file paths, combined table
    without commented rows)`. Single-file mode keeps `(path, table)`.
  - **Loading** (no `-o`/`-a`) still reads a single file, which works for a
    per-setup file through `read_metadata`.

### Command line: `pipnick/scripts/metadata.py`

- New `--split` flag: "Write one metadata file per instrument setup (read
  speed, binning, and readout window), named
  `<base>_<xbin>x<ybin>_<readmode>_<window>`."
- `-t/--raw_table` is the base name in this mode.
- The existence check adapts: in split mode, if any per-setup files exist for
  the base name, `-o` or `-a` is required.

### Reduction: `pipnick/pipelines/reduction.py`

Nothing changes in the code. `pipnick_reduce -t <per-setup file>` reads it
through `build_metadata` → `read_metadata` with the columns restored, and the
setup guard sees one setup. The `rdxdir` update keeps the compact layout.
Update the `reduce_all` docstring to say that a per-setup file can be given.

### Docs and tests

- **Docs:**
  - add a "Per-setup files" section to the `pipnick/metadata.py` module
    docstring;
  - add a release-notes entry to `docs/releases/0.1.0dev.rst`.
- **Tests** (`pipnick/tests/test_metadata.py`, synthetic FITS files, every
  assertion with a message):
  - split mode writes the expected file names, including `full`, `win00` and
    `win01`, and the same window in two setups shares its name;
  - the compact columns and `meta['setup']` are as described;
  - reading back restores the columns, their order, and their types;
  - `expand=False` returns the compact layout;
  - commented rows are preserved per file;
  - a setup with only rejected frames still gets a file;
  - append in split mode adds to the right file, creates a new file with the
    next label and window name for a new setup, and is a no-op when repeated;
  - overwrite warns about stale files;
  - `update_metadata_meta` keeps the compact layout.
- **Test in `pipnick/tests/test_reduction.py`:** a per-setup file passes the
  setup guard.
- **Code style:** follow `AGENTS.md`. Use explicit condition checks, no nested
  functions, and the packed wrapping style. Don't use black.

## Critical files

- `pipnick/metadata.py`: the new functions and the `build_metadata` refactor.
- `pipnick/scripts/metadata.py`: the `--split` flag and the existence check.
- `pipnick/pipelines/reduction.py`: docstring only.
- `pipnick/tests/test_metadata.py`, `pipnick/tests/test_reduction.py`
- `docs/releases/0.1.0dev.rst`
- **Reused:** `assign_setups`, `build_camera_datatable`, `write_metadata`,
  `read_metadata`, `summarize`, `default_table_name` and `_setup_label`, all
  in `pipnick/metadata.py`.

## Verification

- `tox -e 3.13-test,codestyle` passes.
- Run `pipnick_metadata --split --classify stats`, writing to the scratchpad
  only and reading only FITS files from the data directories:
  - **Andor n3:** six files, `n3_1x1_slow_full`, `n3_1x1_fast_full`,
    `n3_2x2_slow_full`, `n3_2x2_fast_full`, `n3_4x4_slow_full` and
    `n3_4x4_fast_full`, with frame counts matching the earlier summary (17/10/
    101/139/10/10 including commented rows).
  - **Andor n2:** `n2_2x2_fast_full` (168 frames) and `n2_2x2_fast_win00` (the
    windowed focus frame).
  - **Old camera n1:** `n1_2x2_slow_full` (61 frames), and
    `n1_2x2_unknown_full` with its single commented row.
- The split files for a night, read back and stacked, match the single-file
  table for the same night: same rows, values and commented rows.
- Append on real data, with scratch symlinks: build split files from part of
  a night, edit one file (comment out a row, change a frame type), link the
  remaining frames, and append. The edits survive, the counts are right, and a
  second append is a no-op.
- `pipnick_reduce -t n1_2x2_slow_full.ecsv` (old camera n1, output to the
  scratchpad) runs end to end as in the previous smoke test. The file keeps
  its compact layout after `rdxdir` is added.

## Progress

### Implementation (done)

- `pipnick/metadata.py`:
  - added `setup_file_keys`, `setup_file_name`, `split_by_setup`,
    `find_setup_files`, `_expand_setup` and `_read_setup_files`;
  - `assign_setups` now assigns `window_name`;
  - `read_metadata` has an `expand` option;
  - `update_metadata_meta` keeps the compact layout;
  - `build_metadata` has a `split` option.
- `pipnick_metadata --split`; `-t` is the base name in this mode.
- **Same rule in the function and the script:** `build_metadata(split=True)`
  raises if per-setup files already exist and neither `overwrite` nor `append`
  is set, as the script does.
- **Small refactor:** `reduce_all`'s setup guard moved, unchanged, into
  `pipnick.pipelines.reduction.check_single_setup`, so it can be tested
  without running a reduction. Its error message now also suggests
  `pipnick_metadata --split`.
- `assign_setups` now adds `window_name` to every setup definition, so
  single-file tables include it in `meta['setups']` too. `test_setups` was
  updated to match.
- Bug fixed during testing: `np.full(n, value, dtype=str)` truncates strings
  to one character. Restored string columns now take their width from the
  value.
- Tests: 61 pass, including 6 new split-mode tests and a per-setup guard test.
  Docs build; release notes updated.

### Verification on real data (scratch outputs only; FITS files only)

- **Andor n3:** `n3_{1x1,2x2,4x4}_{slow,fast}_full`, with 17/10/101/139/10/10
  frames.
- **Andor n2:** `n2_2x2_fast_full` (168) and `n2_2x2_fast_win00` (1).
- **C2 n1:** `n1_2x2_slow_full` (61) and `n1_2x2_unknown_full` (1, commented).
- For all three nights, the split files read back and combined match the
  single-file tables exactly: rows, values and commented rows.
- **Append** (Andor n2 via symlinks, first 100 frames then the rest):
  - 169 unique rows, 22 commented;
  - the user's comment-out and frame-type edit are preserved, and `n2104` is
    auto-commented;
  - a second append leaves the files byte-identical.
- `pipnick_reduce -t n1_2x2_slow_full.ecsv` produced 31 products (bias, four
  flats, 26 science), and the file kept its compact layout after `rdxdir` was
  added.

### Addenda (done)

The user added two behaviors for `--split`. Decisions made with the user:

1. **Pointing and focus frames are commented out by default** in per-setup
   files. The note reads "`<type>` frames are excluded by default". On append,
   this applies only to new frames, so a frame the user restored stays
   restored. Single-file tables are unchanged.
2. **A `comb_id` column** identifies the frames to combine. It is integer, in
   per-setup files only, and placed after `frametype`.
   - **Defaults, within a setup:**
     - all biases share one group;
     - dome flats and sky flats are grouped separately, by filter;
     - darks are grouped by exposure time, rounded to 0.1 s, because the Andor
       `XPOSURE` values jitter by a few ms (e.g., 120.004–120.006 s);
     - all other frames get their own group, including science, focus,
       pointing, test and unknown frames.
   - **Numbering:** IDs are unique across the night, not per file, and
     assigned in file order.
   - **Commented-out frames** also get IDs, so un-commenting needs no other
     edit.
   - **Append:** existing IDs, including user edits, are kept. A new frame
     joins the existing group with the same default key if there is exactly
     one such group. If the user split that group, the new frame gets a new
     group.
   - `pipnick_reduce` does not yet use `comb_id`. That is deferred to a later
     step; the defaults reproduce `reduce_all`'s current behavior.
- **Functions:** `assign_comb_ids` and `_comb_key` add the groups;
  `_exclude_by_default` comments out pointing and focus frames. Both run in
  `build_metadata` when `split=True`, on the new rows only.
- **Tests:** 63 pass. The new tests are `test_comb_ids` and
  `test_comb_ids_append`.
- **Real data, scratch output:**
  - **Andor n4:** bias groups per read mode; sky and dome flats grouped by
    filter (B/V/R/I); darks grouped at 30, 120, 480 and 960 s, with three
    frames each; 13 pointing and focus frames commented out.
  - **C2 n1:** 11 biases in one group; sky flats in groups per filter; 26
    science frames each in their own group; 8 pointing and focus frames
    commented out.
