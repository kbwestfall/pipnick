# Nickel Andor camera: recommended FITS header changes

This document lists recommended changes to the FITS headers written by the
new Nickel Andor science camera. It combines two sources:

- the development team's notes in `Header_metadata.pdf` (the "PDF");
- what the `pipnick` reduction software needed to reduce the 2026
  commissioning data.

Every recommendation was checked against the headers the cameras actually
write:

- **Andor:** all 739 raw frames from the 2026-08-17 run (nights n1–n4);
- **C2:** all 246 raw frames from the 2024-06-24 run (nights n1–n2).

`pipnick` is the new data-reduction software for the Nickel, currently at
https://github.com/kbwestfall/pipnick; it will eventually move to the
UCObservatories organization.

The companion document [`andor_default_wcs.md`](andor_default_wcs.md) gives
the default WCS keywords to insert in each raw frame (item II.1).

## Overview

The recommendations fall into three parts, which have different purposes and
can be scheduled separately:

| Part | Purpose | Who benefits |
|---|---|---|
| [I. Consistency with C2](#i-consistency-with-c2) | Write what the C2 headers provided, with the same names and meanings | Existing user pipelines (Foley, Taggart, Filippenko et al.), unchanged |
| [II. Corrections](#ii-corrections) | Fix Andor values that are wrong or meaningless today | Everyone |
| [III. Improvements beyond C2](#iii-improvements-beyond-c2) | Information C2 never provided, but that makes the reduction more reliable and less camera-specific | New and updated pipelines, including `pipnick` |

Part I is the minimum for existing pipelines to keep working. Part III is
optional, but its first items address the biggest sources of misreduced data
today.

| # | Item | Part | Priority |
|---|---|---|---|
| I.1 | `EXPTIME`: commanded exposure time, 0 for biases | I | **High** |
| I.2 | `OBSTYPE`: `OBJECT` instead of `'?'` | I | **High** |
| I.3 | `DATASEC` (and `COVER`) | I | **High** |
| I.4 | `FILTNAM` (and `FILTRAW`, `FILTORD`) | I | **High** |
| I.5 | `INSTRUME`, `VERSION`, `DSENSOR`, `CAMERAID` | I | **High** |
| I.6 | `HA`, `EQUINOXU` | I | **High** |
| I.7 | `DATE` | I | Medium |
| I.8 | Other C2 keywords: `TEMPDET`, `PROGRAM` | I | Low |
| II.1 | Sky WCS per binning ([`andor_default_wcs.md`](andor_default_wcs.md)); detector WCS follows the binning and the window | II | **High** |
| II.2 | Define what `RA`/`DEC` are | II | Medium |
| II.3 | Placeholder values (`'?'`, blank) | II | Medium |
| III.1 | Explicit frame type chosen by the observer (`IMAGETYP`) | III | **High** |
| III.2 | `BIASSEC`, `CCDSEC` | III | **High** |
| III.3 | `GAIN`, `RDNOISE`, `READMODE`, `SATURATE` | III | **High** |
| III.4 | `SHUTTER` and lamp state | III | Medium |
| III.5 | `CCDSUM` | III | Medium |
| III.6 | `DATE-OBS`, `MJD-OBS`, `TIMESYS` | III | Medium |
| III.7 | `DETECTOR`, cooling set point | III | Low |
| III.8 | Reset `OBJECT` when the frame type changes | III | Low |
| III.9 | Clean-ups: units, the `EXPOSURE` hack, `RADESYS = 'ICRS'` | III | Low |

Also see the [keyword change log](#change-log-c2-to-andor) and the
[decisions for the team](#decisions-for-the-team).

## I. Consistency with C2

The PDF lists the keywords that current users and their pipelines require.
The table below compares each one between C2 and the Andor camera today.

| Keyword | C2 | Andor today | Action |
|---|---|---|---|
| `OBJECT` | Free text | Free text | None |
| `FILTNAM` | Filter name | **Missing**; the name is in `FILTER` | I.4 |
| `INSTRUME` | `'Nickel Direct Camera'` | `'Direct'` | I.5 |
| `TELESCOP` | **Missing** | `'Nickel'` | None; Andor is ahead |
| `VERSION` | `'nickel_direct'` | **Missing** | I.5 |
| `DATASEC` | `'[1:1024,1:1024]'` | **Missing** | I.3 |
| `NAXIS1/2` | Present | Present | None |
| `DATE-BEG` | Exposure start | Exposure start | None; the PDF agrees on `DATE-BEG` over `DATE-STA` |
| `DATE-END` | Exposure end | Exposure end | None |
| `DATE` | UT of the readout | **Missing** | I.7 |
| `EXPTIME` | Commanded time, "not counting shutter error"; 0 for biases | **Missing**; `XPOSURE`, `TELAPSE`, `EXPOSURE` hold the measured time | I.1 |
| `OBSTYPE` | `'DARK'` or `'OBJECT'` | `'DARK'` or `'?'` | I.2 |
| `RA`, `DEC` | Present | Present | Defined in II.2 |
| `HA` | Present | **Missing** | I.6 |
| `AIRMASS` | At the start of the observation | At the start of the observation | None |
| `EQUINOXU` | Equinox of the POCO coordinates | **Missing** | I.6 |
| `EQUINOXS`, `EQUINOX` | Present | Present | None |

### I.1. `EXPTIME`

- C2's `EXPTIME` is the **commanded** exposure time ("not counting shutter
  error"), so it is exactly 0 for biases. Legacy pipelines use it to find
  biases and to identify acquisition images (PDF).
- The Andor camera writes only the measured time. Its biases measure about
  0.08 s, so even a measured-time `EXPTIME` would break the legacy bias
  identification.

**Recommendation.** Write `EXPTIME` as the commanded exposure time in
seconds, **0 for biases**, as C2 did. Keep `XPOSURE` and `TELAPSE` as the
measured values. This resolves the PDF's conflict: the FITS-standard
keywords stay authoritative for the measured time, and `EXPTIME` keeps its
C2 meaning.

### I.2. `OBSTYPE`

- C2 writes only `'DARK'` or `'OBJECT'`, and not consistently. All 11 biases
  on 2024 n1 are `'DARK'`, but all 11 on n2 are `'OBJECT'`. Legacy
  pipelines therefore cannot rely on it alone; the PDF's question of how C2
  sets it remains open.
- The Andor camera writes `'DARK'` when the shutter is closed (as C2 did on
  n1) and `'?'` otherwise.

**Recommendation.** Write `'OBJECT'` instead of `'?'`. Then the Andor values
are a consistent version of C2's: `'DARK'` for a closed shutter, `'OBJECT'`
for an open one.

A **meaningful** frame type is beyond what C2 provided; see III.1.

### I.3. `DATASEC` (and `COVER`)

- C2 writes `DATASEC = '[1:1024,1:1024]'` and `COVER = 32`, the number of
  overscan columns to the right of the data.
- The Andor header's own comments say `DATASEC` is missing.
- Unlike C2, the Andor camera also has a **prescan**, so its data do not
  start in column 1.

**Recommendation.**

- Write `DATASEC` in the binned pixels of the file as written, computed from
  `WINDOW` and the binning. The values are in the table below.
- Write `COVER` only if a user pipeline is found that uses it; it is not in
  the PDF's required list. Note that C2's overscan is clean, but Andor's
  first few overscan columns trail charge. A pipeline that takes the last
  `COVER` columns would include them, which is why III.2 adds `BIASSEC`.

**Measured geometry.** `pipnick` measured this from the 2026 biases and dome
flats, and confirmed it at 1×1, 2×2 and 4×4 and in a window. In unbinned
1-based columns:

- 1–19 are prescan;
- 20–2067 are data;
- 2068–2088 are overscan.

All 2048 rows are data. Binned columns that straddle a boundary are
excluded. The first ~6 overscan columns trail charge, so the bias level is
measured from columns 2075–2088 (see III.2). The resulting values, which can
serve as test cases:

| Binning | Window (`WINDOW`) | Shape | `DATASEC` | `BIASSEC` (III.2) |
|---|---|---|---|---|
| 1×1 | full | 2088×2048 | `[20:2067,1:2048]` | `[2075:2088,1:2048]` |
| 2×2 | full | 1044×1024 | `[11:1033,1:1024]` | `[1038:1044,1:1024]` |
| 4×4 | full | 522×512 | `[6:516,1:512]` | `[520:522,1:512]` |
| 2×2 | `960 1205 840 1069` | 123×115 | `[1:123,1:115]` | none (omit) |

### I.4. `FILTNAM`

- C2 writes `FILTNAM`, plus the wheel position as `FILTRAW` and `FILTORD`.
- The Andor camera writes `FILTER` instead.
- The PDF lists `FILTNAM` as required, and notes that some users would
  prefer `FILTER`.

**Recommendation.**

- Write `FILTNAM` and `FILTER` with the same value, so that both groups of
  users are served.
- Add `FILTRAW` and `FILTORD` if the wheel position is available.
- **Beware:** consider avoiding quote characters in filter names (e.g. `i'`).
  They are legal in FITS (written as `i''`), but they break naive parsers and
  shell quoting.

### I.5. Instrument identification

- C2 identifies itself with `INSTRUME = 'Nickel Direct Camera'`,
  `VERSION = 'nickel_direct'`, `DSENSOR = 'Loral 2Kx2K'` and `CAMERAID = 2`.
- The Andor camera writes `INSTRUME = 'Direct'`, no `VERSION`, `DSENSOR` or
  `CAMERAID`, and identifies itself only in the free-form `CAMID`.
- The PDF says pipelines use `INSTRUME`, `TELESCOP` and `VERSION` to choose
  their reduction path, expects a new value for the Andor camera, and wants
  all three kept.

**Recommendation.**

- `INSTRUME`: a stable value that distinguishes the Andor camera from C2.
- `VERSION`: a new value for the Andor data-acquisition system.
- `DSENSOR`: e.g. `'Andor DZ936_BV'`.
- `CAMERAID`: a new number.

These let pipelines select the camera from the same keywords they use for C2.

### I.6. `HA` and `EQUINOXU`

- C2 writes `HA`, which the PDF says pipelines use to sort sky flats. The
  Andor camera does not.
- C2 writes `EQUINOXU`, the equinox of the POCO coordinates; its comment
  reads "epoch for POCO position is current date" (e.g., `2024.48`). The
  PDF says one team relies on it exclusively. The Andor camera does not write
  it.

**Recommendation.**

- Add `HA`.
- Add `EQUINOXU`, giving the equinox of `RA`/`DEC`.
  - The 2026 data show that the Andor `RA`/`DEC` are J2000: they agree with
    the astrometric solutions to within 30″, whereas coordinates of date
    would be off by several arcminutes after 26 years of precession (see
    [`andor_default_wcs.md`](andor_default_wcs.md)).
  - So `EQUINOXU = 2000.0`, unless the telescope-control system says
    otherwise (II.2).

### I.7. `DATE`

C2 writes `DATE` as the UT of the readout; the PDF lists it as required.

**Recommendation.** Add `DATE`, the time the file was written. This is the
FITS-standard meaning, and it is effectively the same as C2's.

### I.8. Other C2 keywords (low priority)

The PDF does not list these as required, but they have direct Andor
equivalents under other names. Add them only if a user pipeline turns out to
need them:

- `TEMPDET`, the detector temperature (Andor: `CURRTEMP`);
- `PROGRAM` (Andor: `PROGNAME`, always blank in 2026).

Some C2 keywords have **no** meaningful Andor equivalent and should *not* be
imitated:

- `GAIN`: C2's is a DCS index, not e⁻/ADU. Do not write a C2-style index; see
  III.3.
- `READ-SPD`: C2's readout-speed code. However, consider defining an
  enumeration for the `AMPCONF` values and adding it; cf. `READMODE` in III.3.

## II. Corrections

These Andor values are wrong or ambiguous today.

### II.1. The detector and sky WCS must follow the binning and the window

Every 2026 frame carries the same detector (`…C`), sky (`…S`) and primary
WCS, whatever its binning or window:

- `CD1_1C = 2.0` and `CRVAL1C = -20.0`;
- `CRPIX1/2 = 522, 513`;
- the 2×2 plate scale.

These values are correct only for full-frame 2×2, and wrong for the 1×1,
4×4 and windowed frames.

**Recommendation.**

- **Sky WCS (primary and `…S`):** use the static keyword sets in
  [`andor_default_wcs.md`](andor_default_wcs.md). There is one set for each
  binning (1×1, 2×2, 4×4), measured from the 2026 data against Gaia.
  - `CRPIX` is the center of the binned frame.
  - `CRVAL` is `RA`/`DEC`.
  - The `CD` matrix is either the measured one, or an orthogonal one with a
    rotation of exactly 180°.

  That document also gives:
  - the uncertainties of the keywords;
  - suggested `CRDER`/`CSYER` values;
  - a comparison with the current default.

  This header WCS is only a first guess. Observers must still solve the
  astrometry themselves (e.g., with `pipnick`) to get accurate coordinates.
- **Windowed frames:** the static sets are valid only for unwindowed
  readouts. For a windowed readout, omit the sky WCS or flag it as invalid,
  unless the header writer can shift `CRPIX` for the window origin.
- **Detector WCS:** `CD1_1C`/`CD2_2C` are the binning, and `CRVAL1C`/`CRVAL2C`
  include the window origin.

**Also check:** the detector WCS (`CRVAL1C = -20`, `CRPIX1C = 0.5`) puts
detector coordinate 0 at the left edge of unbinned column 21, but the data
begin at column 20. The detector WCS should agree with `DATASEC`.

### II.2. Define what `RA`/`DEC` are

**Current state.**

- The comments say `Right Ascension if POCOT=0`, while `RA_T` is
  `Right Ascension if POCOT!=0`, and another comment says the coordinates
  "are intended for engineering, not science". None of this says what the
  value is.
- In the 2026 data, `RA`/`DEC` agree with the `pipnick` astrometric
  solutions to within 30″ (median 22″), for both `POCOT` values; see
  [`andor_default_wcs.md`](andor_default_wcs.md).
- `RA_D`/`RA_T` are the target before dither offsets, which put them 100″ to
  135″ off for dithered frames.

**Recommendation.**

- Make `RA`/`DEC` the telescope's actual pointing, as on C2, with an
  unambiguous comment.
- State their equinox in `EQUINOXU` (I.6). The data indicate J2000.

### II.3. Update keyword values that are currently `?` or remove them

Some header keywords carry no information. In the 2026 data, these are:

- `COOLING = '?'` in every frame;
- `OBSTYPE = '?'` in every open-shutter frame (see I.2);
- `PROGNAME` blank in every frame.

These should be updated with real values or removed.

## III. Improvements beyond C2

C2 never provided these, so existing pipelines don't need them. They make the
reduction more reliable and less dependent on camera-specific rules; today,
`pipnick` has to infer each of them.

### III.1. Explicit frame type chosen by the observer

Neither camera records what a frame *is*, so it has to be inferred from the
shutter state, the exposure time and the free-text `OBJECT`, which is often
stale. This is the largest source of misclassified frames. In the 2026 data:

- biases labelled 'test';
- open-shutter 'test bias' frames;
- flats keeping an earlier 'twilight flat' label.

**Recommendation.**

- The observing software sets `IMAGETYP` from the observer's selection,
  using a fixed vocabulary:

  ```
  bias, dark, domeflat, skyflat, object, focus, pointing, test
  ```
- Leave `OBSTYPE` with its C2 meaning (I.2), so that legacy pipelines are
  unaffected.

### III.2. `BIASSEC` and `CCDSEC`

- **`BIASSEC`:** the overscan to use for the bias level. C2 gave only
  `COVER`. For Andor, it should exclude the columns that trail charge; see
  the table in I.3. Omit it when a window contains no overscan.
- **`CCDSEC`:** the unbinned detector pixels covered by `DATASEC`, so that
  binned and windowed frames map back to the detector. This matters for
  bad-pixel masks and for matching calibrations between setups. The IRAF
  `LTV1/2` and `LTM1_1/LTM2_2` are an acceptable alternative.
- Optional: `PRESEC` for the prescan, and `DETSIZE = '[1:2088,1:2048]'`.

### III.3. Gain, read noise, readout mode and saturation

**Current state.** The readout mode is recorded only as `AMPCONF`
(`'1.0MHz'` or `'0.05MHz'`) and numeric IDs (`PREAMPID`, `AMPID`, `HSSID`).
`pipnick` has to look up the gain and read noise.

**Recommendation.**

- Add `GAIN` [e⁻/ADU] and `RDNOISE` [e⁻] for the mode actually used, from a
  table maintained with the camera software. Andor's performance-report
  values:

  | `AMPCONF` | Gain | Read noise |
  |---|---|---|
  | 1.0 MHz | 0.91 e⁻/ADU | 6.09 e⁻ |
  | 0.05 MHz | 0.92 e⁻/ADU | 2.70 e⁻ |

  C2's `GAIN` has a different meaning (a DCS index); see the
  [change log](#change-log-c2-to-andor).
- Add `READMODE = 'fast'/'slow'` alongside `AMPCONF`; cf. I.8.
- Add `SATURATE` [ADU]. The nominal saturation is the standard `2**16-1`
  (65535) ADU, although the observed maxima are about 65240 ADU at 1.0 MHz
  and 65140 ADU at 0.05 MHz.

### III.4. Shutter and lamp state

- **`SHUTTER`:** add it, with a fixed vocabulary (e.g. `open`, `closed`).
  Today the shutter state is available only implicitly through
  `OBSTYPE = 'DARK'`.
- **Lamps:** `LAMPS = 'Off'` in every frame, including dome flats (e.g.
  n2100, `OBJECT = 'Domeflat B'`), so it evidently reports only the Hg/Ne arc
  lamps.
  - If the dome-flat lamps are under software control, record their state.
  - Otherwise, the `LAMPS` comment should say it covers only the arc lamps.

### III.5. `CCDSUM`

Both cameras record the binning only in private UCAM/xvideo keywords:

- `CDELT1U/2U`, whose sign differs between C2 (−2) and Andor (+2);
- `CBIN/RBIN` on Andor, which are log₂ values.

**Recommendation.** Add `CCDSUM = '2 2'` (the IRAF convention: column binning
first, then row binning).

### III.6. Standard time keywords

Add these, which most general tools read:

- `DATE-OBS`, with the same value as `DATE-BEG`;
- `MJD-OBS` (or `MJD-BEG`, plus `MJD-END`);
- `TIMESYS = 'UTC'`.

### III.7. `DETECTOR` and cooling set point

- **`DETECTOR`:** e.g. `'Andor DZ936_BV SN 31382'`; the standard
  counterpart of `DSENSOR` (I.5).
- **Cooling set point** (e.g. `CCDTSET`) next to `CURRTEMP`. `pipnick` rejects
  warm frames (e.g. n1017–1018), and the set point makes "warm" well defined.

### III.8. Reset `OBJECT` when the frame type changes

Stale `OBJECT` labels mislead both human readers and any pipeline that reads
them. Example: n4046–4055, labelled 'twilight flat' at 384–727 ADU. With III.1
in place, this matters less.

### III.9. Clean-ups

- **Units in comments:** units go first, in brackets (PDF). This is done for
  most keywords; the main one remaining is `VSSPEED` (presumably µs per row).
- **`EXPOSURE`:** its comment calls it `Hack: same as XPOSURE`. Remove it
  once `EXPTIME` (I.1) exists.
- **Private keywords:** these are now labelled `UCAM private`,
  `XVIDEO private` or `Engr private` (PDF). Keep them, since xvideo needs
  `CRVALiU`/`CDELTiU`.
- **`[XY]PIXSIZE`, `[XY]DETSIZE`:** keep them. The PDF notes that they aren't
  Andor-specific.
- **`RADESYS = 'ICRS'`** (with no `EQUINOX`) for the WCS, since the
  astrometric references (Gaia, astrometry.net) are ICRS. `FK5` J2000
  differs by less than 0.1″, so this is low priority. C2 wrote the deprecated
  `RADECSYS`; Andor's `RADESYS` is already correct.

## Change log: C2 to Andor

The keywords that pipelines written for C2 use, and how each changes for the
Andor camera once the recommendations above are adopted. Keywords that are
new on the Andor camera, with no C2 counterpart, are listed separately below.

| C2 keyword | C2 value or meaning | Andor keyword(s) | Andor value or meaning | Item |
|---|---|---|---|---|
| `OBJECT` | Free text | `OBJECT` | Unchanged; reset when the frame type changes | III.8 |
| `OBSTYPE` | `'DARK'` or `'OBJECT'` | `OBSTYPE` | Unchanged | I.2 |
| `EXPTIME` | Commanded time [s]; 0 for biases | `EXPTIME` | Unchanged. The measured time is in `XPOSURE`, `TELAPSE` | I.1 |
| `FILTNAM` | Filter name | `FILTNAM`, `FILTER` | Unchanged; also written as `FILTER` | I.4 |
| `FILTRAW`, `FILTORD` | Filter-wheel position | `FILTRAW`, `FILTORD` | Unchanged, if available | I.4 |
| `INSTRUME` | `'Nickel Direct Camera'` | `INSTRUME` | **New value**, identifying the Andor camera | I.5 |
| `VERSION` | `'nickel_direct'` | `VERSION` | **New value** | I.5 |
| `DSENSOR` | `'Loral 2Kx2K'` | `DSENSOR`, `DETECTOR` | **New value**, e.g. `'Andor DZ936_BV'` | I.5, III.7 |
| `CAMERAID` | `2` | `CAMERAID` | **New value** | I.5 |
| `DATASEC` | `'[1:1024,1:1024]'` | `DATASEC` | **Shifted by the prescan**; depends on the binning and window (e.g. `'[11:1033,1:1024]'` at 2×2) | I.3 |
| `COVER` | `32` overscan columns | `BIASSEC` (`COVER` only if needed) | The overscan section to use, excluding columns that trail charge | I.3, III.2 |
| `DATE` | UT of the readout | `DATE` | Time the file was written | I.7 |
| `DATE-BEG`, `DATE-END` | Exposure start and end | `DATE-BEG`, `DATE-END` | Unchanged | — |
| `RA`, `DEC` | Pointing | `RA`, `DEC` | Actual pointing, J2000 | II.2 |
| `EQUINOXU` | Equinox of `RA`/`DEC` (current date) | `EQUINOXU` | Equinox of `RA`/`DEC` (2000.0) | I.6 |
| `HA` | Hour angle | `HA` | Unchanged | I.6 |
| `AIRMASS` | At the start of the observation | `AIRMASS` | Unchanged | — |
| `GAIN` | **DCS gain index** | `GAIN` | **e⁻/ADU**, for the readout mode used | III.3 |
| `READ-SPD` | DCS readout-speed code (80) | `AMPCONF`, `READMODE` | `'1.0MHz'`/`'0.05MHz'`; `'fast'`/`'slow'` | III.3 |
| `TEMPDET` | Detector temperature [°C] | `CURRTEMP` (`TEMPDET` only if needed) | Unchanged meaning | I.8 |
| `PROGRAM` | Program name | `PROGNAME` (`PROGRAM` only if needed) | Unchanged meaning | I.8 |
| `CDELT1U/2U` | −2 at 2×2 (private) | `CDELT1U/2U` | **+2** at 2×2 (private; use `CCDSUM`) | III.5 |
| `CRVAL1U/2U` | 2048 (private) | `CRVAL1U/2U` | 0 (private) | — |
| `RADECSYS` | `'FK5'` (deprecated keyword) | `RADESYS` | `'FK5'` (or `'ICRS'`) | III.9 |
| `CRPIX1/2` | 512, 512 | `CRPIX1/2` | Center of the binned frame, e.g. 522.5, 512.5 at 2×2 | II.1 |
| `CD1_1` … `CD2_2` | 0.3698″/pix at 2×2, rotation 182.2° | `CD1_1` … `CD2_2` | 0.3331″/pix at 2×2, rotation 180.55° (or 180°), per binning | II.1 |
| `CD1_1C`, `CRVAL1C` (detector WCS) | −2, 2048 | `CD1_1C`, `CRVAL1C` | +binning; origin at the prescan, e.g. −20 | II.1 |

**New on the Andor camera, with no C2 counterpart:**

- already present: `TELESCOP`, `XPOSURE`, `TELAPSE`, `WINDOW`, `AMPCONF`;
- recommended: `IMAGETYP`, `BIASSEC`, `CCDSEC`, `RDNOISE`, `READMODE`,
  `SATURATE`, `SHUTTER`, `CCDSUM`, `DATE-OBS`, `MJD-OBS`, `TIMESYS`,
  `DETECTOR`, and a cooling set point.

## Decisions for the team

1. **`EXPTIME`.** The PDF says it "will not be added", but it is also listed
   as required and as a TODO. Recommendation: the commanded time, as on C2
   (I.1).
2. **`FILTER` vs `FILTNAM`.** Recommendation: write both (I.4).
3. **`INSTRUME`, `VERSION`, `CAMERAID` values for the Andor camera** (I.5).
4. **Measured or orthogonal CD matrix** for the default WCS (II.1;
   [`andor_default_wcs.md`](andor_default_wcs.md)). Either is adequate as a
   first guess.
5. **Whether to adopt Part III now**, and in particular III.1 (`IMAGETYP`),
   which needs a change in the observing interface as well as in the header
   writer.

## Already done

These are done in the 2026 headers, since the PDF's v2 example:

- Units in brackets at the start of most comments.
- Private UCAM/xvideo/engineering keywords labelled as such.
- `XPOSURE` and `TELAPSE` added; `TIMEUNIT = 's'`.
- `AIRMASS` added, at the start of the exposure, as on C2.
- `EQUINOX` and `RADESYS` present, as part of a new default primary WCS. The
  bogus `EQUINOXC = 0.0` from v2 is gone.
- A default sky WCS (primary and alternate `S`) added, "to be replaced by
  astrometry". This has the same structure as C2. Keeping `S` as an
  untouched copy is good practice, since astrometric pipelines overwrite the
  primary. The updated values are in
  [`andor_default_wcs.md`](andor_default_wcs.md) (II.1).
- `TELESCOP = 'Nickel'`, which C2 never had.
- `WINDOW` comment clarified (`Andor window: hbeg hend vbeg vend`; 1-based,
  unbinned).
