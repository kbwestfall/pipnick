# Nickel Andor camera: static default WCS

This document provides the default (pre-astrometry) sky WCS keywords that
the Andor header writer should insert in each raw frame, one set for each
binning, measured from the 2026-08-17 commissioning data. It accompanies
[`andor_header_recommendations.md`](andor_header_recommendations.md)
(item II.1).

## Recommended keywords

Each block below is exactly what the header writer should insert as the
primary WCS of an **unwindowed** frame. The header writer chooses the block
by binning. Adopt **either** the measured or the orthogonal CD matrix
(see [Alternative: an orthogonal WCS](#alternative-an-orthogonal-wcs-rotation-of-exactly-180)).

- **`CRVAL1/2`:** these depend on the frame; they are the `RA`/`DEC` keywords
  converted to degrees. The blocks show them for an example pointing,
  `RA = '22:37:04.00'` and `DEC = '34:24:59.0'` (n2117, NGC 7331).
- **Everything else:** fixed for a given binning.
- **The alternate `…S` WCS:** use the same values, with the `S` suffix on
  each keyword (`CRPIX1S`, `CD1_1S`, `RADESYSS`, etc.).
- **Construction:** the blocks were built with `astropy.io.fits`, and read
  back with `astropy.wcs`. This confirmed the scale (0.16655″ per unbinned
  pixel), the position angle (180.55° measured, 180.00° orthogonal) and
  that the reference pixel falls on `CRVAL`.

### Measured CD matrix (rotation 180.55°)

1×1 (2088 × 2048):

```
WCSNAME = 'Celestial coordinates' / name of WCS
CTYPE1  = 'RA---TAN'           / type of WCS axis i
CTYPE2  = 'DEC--TAN'           / type of WCS axis i
CUNIT1  = 'deg     '           / unit of WCS axis i
CUNIT2  = 'deg     '           / unit of WCS axis i
CRPIX1  =               1044.5 / [FITSpix] reference pixel: frame center
CRPIX2  =               1024.5 / [FITSpix] reference pixel: frame center
CRVAL1  =           339.266667 / [deg] RA keyword at reference pixel
CRVAL2  =            34.416389 / [deg] DEC keyword at reference pixel
CD1_1   =          -4.6262E-05 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD1_2   =            -4.44E-07 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD2_1   =            4.386E-07 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD2_2   =          -4.6261E-05 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
RADESYS = 'FK5     '           / reference system for celestial coordinates
EQUINOX =               2000.0 / [a] date of WCS celestial reference frame
CRDER1  =              0.00025 / [deg] random error along WCS axis i
CRDER2  =              0.00025 / [deg] random error along WCS axis i
CSYER1  =               0.0083 / [deg] systematic error along WCS axis i
CSYER2  =               0.0083 / [deg] systematic error along WCS axis i
```

2×2 (1044 × 1024):

```
WCSNAME = 'Celestial coordinates' / name of WCS
CTYPE1  = 'RA---TAN'           / type of WCS axis i
CTYPE2  = 'DEC--TAN'           / type of WCS axis i
CUNIT1  = 'deg     '           / unit of WCS axis i
CUNIT2  = 'deg     '           / unit of WCS axis i
CRPIX1  =                522.5 / [FITSpix] reference pixel: frame center
CRPIX2  =                512.5 / [FITSpix] reference pixel: frame center
CRVAL1  =           339.266667 / [deg] RA keyword at reference pixel
CRVAL2  =            34.416389 / [deg] DEC keyword at reference pixel
CD1_1   =          -9.2523E-05 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD1_2   =           -8.879E-07 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD2_1   =            8.773E-07 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD2_2   =          -9.2523E-05 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
RADESYS = 'FK5     '           / reference system for celestial coordinates
EQUINOX =               2000.0 / [a] date of WCS celestial reference frame
CRDER1  =              0.00025 / [deg] random error along WCS axis i
CRDER2  =              0.00025 / [deg] random error along WCS axis i
CSYER1  =               0.0083 / [deg] systematic error along WCS axis i
CSYER2  =               0.0083 / [deg] systematic error along WCS axis i
```

4×4 (522 × 512):

```
WCSNAME = 'Celestial coordinates' / name of WCS
CTYPE1  = 'RA---TAN'           / type of WCS axis i
CTYPE2  = 'DEC--TAN'           / type of WCS axis i
CUNIT1  = 'deg     '           / unit of WCS axis i
CUNIT2  = 'deg     '           / unit of WCS axis i
CRPIX1  =                261.5 / [FITSpix] reference pixel: frame center
CRPIX2  =                256.5 / [FITSpix] reference pixel: frame center
CRVAL1  =           339.266667 / [deg] RA keyword at reference pixel
CRVAL2  =            34.416389 / [deg] DEC keyword at reference pixel
CD1_1   =         -0.000185047 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD1_2   =           -1.776E-06 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD2_1   =            1.755E-06 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD2_2   =         -0.000185045 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
RADESYS = 'FK5     '           / reference system for celestial coordinates
EQUINOX =               2000.0 / [a] date of WCS celestial reference frame
CRDER1  =              0.00025 / [deg] random error along WCS axis i
CRDER2  =              0.00025 / [deg] random error along WCS axis i
CSYER1  =               0.0083 / [deg] systematic error along WCS axis i
CSYER2  =               0.0083 / [deg] systematic error along WCS axis i
```

### Orthogonal CD matrix (rotation exactly 180°)

1×1 (2088 × 2048):

```
WCSNAME = 'Celestial coordinates' / name of WCS
CTYPE1  = 'RA---TAN'           / type of WCS axis i
CTYPE2  = 'DEC--TAN'           / type of WCS axis i
CUNIT1  = 'deg     '           / unit of WCS axis i
CUNIT2  = 'deg     '           / unit of WCS axis i
CRPIX1  =               1044.5 / [FITSpix] reference pixel: frame center
CRPIX2  =               1024.5 / [FITSpix] reference pixel: frame center
CRVAL1  =           339.266667 / [deg] RA keyword at reference pixel
CRVAL2  =            34.416389 / [deg] DEC keyword at reference pixel
CD1_1   =          -4.6264E-05 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD1_2   =                  0.0 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD2_1   =                  0.0 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD2_2   =          -4.6264E-05 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
RADESYS = 'FK5     '           / reference system for celestial coordinates
EQUINOX =               2000.0 / [a] date of WCS celestial reference frame
CRDER1  =              0.00072 / [deg] random error along WCS axis i
CRDER2  =              0.00072 / [deg] random error along WCS axis i
CSYER1  =               0.0083 / [deg] systematic error along WCS axis i
CSYER2  =               0.0083 / [deg] systematic error along WCS axis i
```

2×2 (1044 × 1024):

```
WCSNAME = 'Celestial coordinates' / name of WCS
CTYPE1  = 'RA---TAN'           / type of WCS axis i
CTYPE2  = 'DEC--TAN'           / type of WCS axis i
CUNIT1  = 'deg     '           / unit of WCS axis i
CUNIT2  = 'deg     '           / unit of WCS axis i
CRPIX1  =                522.5 / [FITSpix] reference pixel: frame center
CRPIX2  =                512.5 / [FITSpix] reference pixel: frame center
CRVAL1  =           339.266667 / [deg] RA keyword at reference pixel
CRVAL2  =            34.416389 / [deg] DEC keyword at reference pixel
CD1_1   =          -9.2527E-05 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD1_2   =                  0.0 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD2_1   =                  0.0 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD2_2   =          -9.2527E-05 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
RADESYS = 'FK5     '           / reference system for celestial coordinates
EQUINOX =               2000.0 / [a] date of WCS celestial reference frame
CRDER1  =              0.00072 / [deg] random error along WCS axis i
CRDER2  =              0.00072 / [deg] random error along WCS axis i
CSYER1  =               0.0083 / [deg] systematic error along WCS axis i
CSYER2  =               0.0083 / [deg] systematic error along WCS axis i
```

4×4 (522 × 512):

```
WCSNAME = 'Celestial coordinates' / name of WCS
CTYPE1  = 'RA---TAN'           / type of WCS axis i
CTYPE2  = 'DEC--TAN'           / type of WCS axis i
CUNIT1  = 'deg     '           / unit of WCS axis i
CUNIT2  = 'deg     '           / unit of WCS axis i
CRPIX1  =                261.5 / [FITSpix] reference pixel: frame center
CRPIX2  =                256.5 / [FITSpix] reference pixel: frame center
CRVAL1  =           339.266667 / [deg] RA keyword at reference pixel
CRVAL2  =            34.416389 / [deg] DEC keyword at reference pixel
CD1_1   =         -0.000185054 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD1_2   =                  0.0 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD2_1   =                  0.0 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
CD2_2   =         -0.000185054 / [deg/FITSpix] matrix i_j FITSax j => WCSax i
RADESYS = 'FK5     '           / reference system for celestial coordinates
EQUINOX =               2000.0 / [a] date of WCS celestial reference frame
CRDER1  =              0.00072 / [deg] random error along WCS axis i
CRDER2  =              0.00072 / [deg] random error along WCS axis i
CSYER1  =               0.0083 / [deg] systematic error along WCS axis i
CSYER2  =               0.0083 / [deg] systematic error along WCS axis i
```

## Assumptions

1. **One fixed keyword set per binning.** The header writer inserts the WCS
   by copying a fixed set of keywords, chosen by binning (1×1, 2×2, 4×4).
   - Windowed readouts are **not** supported: their WCS would be wrong.
   - Only 2×2 frames were taken on sky in 2026-08-17. The 2×2 set is
     measured; the 1×1 and 4×4 sets are derived from it.
2. **`CRPIX1/2` is the center of the binned, unwindowed frame**, including
   the prescan and overscan columns. In FITS (1-based) pixels, that is
   `(NAXIS1 + 1)/2, (NAXIS2 + 1)/2`.
3. **`CRVAL1/2` equals the `RA`/`DEC` header keywords**, in degrees. The
   offset between the pointing and the actual field center is ignored when
   deriving the CD matrix. It is quantified below only as the dominant
   uncertainty of `CRVAL`.

## Data and method

- **Frames:** all 53 science frames from 2026-08-17 with a `pipnick`
  astrometric solution. All are 2×2 full-frame:
  - n2: NGC 7331, 15 frames, 1.0 MHz;
  - n3: SA 111 (24 frames) and M11 (14 frames), both readout speeds.
  - Nights n1 and n4 have not been reduced.
- **Measurement:** the CD matrix of each frame was measured independently of
  `pipnick`'s solution, which for most frames is a bootstrap from another
  frame.
  - Detected sources were taken from `pipnick`'s photometry tables, using only
    sources without flags.
  - They were matched to Gaia DR3 (ICRS), propagated to the epoch of the
    observation, with mutual nearest neighbours within 1″. Each frame has
    28–359 matches.
  - The Gaia positions were projected onto the tangent plane at the
    **center of the raw (untrimmed) frame**.
  - A linear model (CD matrix plus a center offset) was fit by least
    squares with 3σ clipping. The tangent point was iterated to the fitted
    center.
- **Distortion check:** each frame was also fit with a cubic model.

## Results (2×2)

**Fit quality.** The linear fits have RMS residuals of 0.03–0.17″, i.e.,
0.1–0.5 binned pixels. The one exception is M11 V C096, at 0.40″.

**No SIP distortion is needed.** In the dense fields (M11, SA 111), the cubic
part of the fit departs from the best linear model by a median of 0.10″. It
is at most 0.21″ anywhere in the frame in all but two frames: 0.54″ (SA 111 I
C045) and 0.90″ (the M11 V C096 outlier). Using
the cubic fit instead of the linear one changes the scale by
0.001% and the rotation by 1″, both negligible.

**Scale and rotation.**

| Sample | Frames | Scale (″/pix) | Position angle of +y (deg, E of N) |
|---|---|---|---|
| All (adopted) | 53 | 0.33310 | 180.550 |
| n2 NGC 7331 | 15 | 0.33314 | 180.366 |
| n3 SA 111 | 24 | 0.33308 | 180.551 |
| n3 M11 | 14 | 0.33308 | 180.592 |
| *Current header default* | — | *0.33307* | *180.300* |

- **Orientation:** the determinant of the CD matrix is positive in every frame.
  With `CD1_1 < 0` and a position angle near 180°:
  - north is toward −y and east toward −x;
  - the image is mirror-reversed relative to the usual north-up,
    east-left orientation.
- **Pixels:** square to 0.01% (per-axis scales agree), and orthogonal to
  0.01°.
- **Scale:** stable to about 0.02% between frames, fields and nights.
- **Rotation:** this is the only quantity that varies significantly.
  - Within a field, the scatter is about 0.01°.
  - Between fields on the same night, it differs by 0.04° (SA 111 vs. M11).
  - Between nights, it differs by 0.19° (n2 vs. n3).
  - No correlation with airmass or time was found within n3 (|r| ≤ 0.3).
  - Whether the night-to-night difference is a change in the camera's
    mounting, flexure, or something else cannot be determined from two
    nights. The rotation should be monitored as more data accumulate.
- **Readout speed:** it makes no difference beyond the field-to-field
  variation.

## Adopted keyword sets

**Values.** The CD matrix is the **element-wise median of all 53 frames**.
The median of the three field medians differs from it by less than a third of
the uncertainties below. The 1×1 and 4×4 matrices are the 2×2 matrix
scaled by 1/2 and 2.

**Common to all three sets:**

```
WCSNAME = 'Celestial coordinates'
CTYPE1  = 'RA---TAN'
CTYPE2  = 'DEC--TAN'
CUNIT1  = 'deg'
CUNIT2  = 'deg'
CRVAL1  = <RA in degrees>
CRVAL2  = <DEC in degrees>
RADESYS = 'FK5'          (or 'ICRS' without EQUINOX; see below)
EQUINOX = 2000.0
```

**Per binning:**

| Keyword | 1×1 (2088×2048) | 2×2 (1044×1024) | 4×4 (522×512) |
|---|---|---|---|
| `CRPIX1` | 1044.5 | 522.5 | 261.5 |
| `CRPIX2` | 1024.5 | 512.5 | 256.5 |
| `CD1_1` | −4.6262E−05 | −9.2523E−05 | −1.85047E−04 |
| `CD1_2` | −4.440E−07 | −8.879E−07 | −1.776E−06 |
| `CD2_1` | 4.386E−07 | 8.773E−07 | 1.755E−06 |
| `CD2_2` | −4.6261E−05 | −9.2523E−05 | −1.85045E−04 |

The CD values are rounded to the precision their uncertainties support.
The reference pixel corresponds to the same detector position, unbinned
column 1044.5 and row 1024.5, in all three binnings.

The same values apply to the alternate `…S` WCS. The detector WCS (`…C`) is a
separate issue; see
[`andor_header_recommendations.md`](andor_header_recommendations.md) (II.1).

## Uncertainties

**CD matrix (2×2; scale the absolute values by 1/2 or 2 for 1×1 and 4×4):**

| Element | Value (deg/pix) | Statistical | Systematic | Adopted uncertainty |
|---|---|---|---|---|
| `CD1_1` | −9.2523E−05 | 5E−09 | 2.1E−08 | **2E−08** (0.02%) |
| `CD1_2` | −8.879E−07 | 1.4E−08 | 1.8E−07 | **1.8E−07** (20%) |
| `CD2_1` | 8.773E−07 | 1.4E−08 | 1.8E−07 | **1.8E−07** (20%) |
| `CD2_2` | −9.2523E−05 | 2.4E−09 | 4.4E−09 | **5E−09** (0.005%) |

- *Statistical* is the error of the median: 1.25 × the robust scatter /
  √53.
- *Systematic* is half the range of the three field medians (n2 NGC 7331,
  n3 SA 111, n3 M11). It dominates.
- **Derived quantities:** the scale is 0.33310 ± 0.00005″/pix (±0.015%), and
  the rotation is 180.55 ± 0.12°.
- The large *fractional* uncertainty of the off-diagonal elements is just
  the rotation uncertainty: those elements are small because the camera is
  nearly aligned with the sky axes.

**Effect on positions.** The error from using the static CD instead of each
frame's own, at the corners of the 2×2 frame (243″ from the center):

| CD matrix | Median | 90th percentile | Max |
|---|---|---|---|
| Adopted | 0.18″ | 0.83″ | 0.86″ |
| Current header default | 1.08″ | 1.26″ | 1.31″ |

The worst cases for the adopted matrix are the n2 frames, whose rotation
differs by 0.19°.

**`CRVAL`, i.e., the pointing.** The field center found by astrometry
differs from `RA`/`DEC`:

- median 22″, 90th percentile 28″, maximum 30″;
- stable to 0.2–1.5″ within a field, but different from field to field. For
  n2 NGC 7331: +23″ in RA × cos Dec and −16″ in Dec; for n3 M11: +20″ and
  +6″.

This is about 25 times larger than the error from the CD matrix, so the
pointing completely dominates the accuracy of the default WCS. If the error
keywords are kept, suitable values for all three binnings would be:

| Keyword | Value | Equivalent | Meaning |
|---|---|---|---|
| `CSYER1/2` | 0.0083 deg | 30″ | Pointing; the current value is 60″, which is conservative |
| `CRDER1/2` | 0.00025 deg | 0.9″ | CD-matrix error at the frame corners; the current value is 0.17″, which is too small |
| `CRDER1/2` (orthogonal CD) | 0.00072 deg | 2.6″ | As above, for the orthogonal matrix |

## Comparison with the current header default

The current default is already close for 2×2 full-frame data. The adopted set
changes:

1. **`CD`:** the rotation goes from 180.30° to 180.55°. The scale is
   essentially unchanged (0.33307 to 0.33310″/pix). The corner errors fall
   from about 1.1″ to a median of 0.2″.
2. **`CRPIX`:** from `(522, 513)` to the exact frame center, `(522.5, 512.5)`.
   This is half a pixel in each axis, i.e., about 0.2″.
3. **Separate sets for 1×1 and 4×4:** the current default is wrong by a
   factor of 2 in scale (and in `CRPIX`) for those binnings.

## Alternative: an orthogonal WCS (rotation of exactly 180°)

The header WCS is only a first guess: observers must solve the astrometry
themselves (e.g., with `pipnick`) to get accurate coordinates. The
off-diagonal terms can therefore be dropped, keeping the adopted scale and
the mirror-reversed orientation (positive determinant).

**Keywords:** `CD1_2 = CD2_1 = 0`, plus the following, with the same
`CRPIX` and `CRVAL` as above:

| Keyword | 1×1 | 2×2 | 4×4 |
|---|---|---|---|
| `CD1_1` | −4.6264E−05 | −9.2527E−05 | −1.85054E−04 |
| `CD2_2` | −4.6264E−05 | −9.2527E−05 | −1.85054E−04 |

**Offset** from each frame's measured CD, evaluated over the whole frame:

| | Orthogonal: median frame | Orthogonal: worst frame | Adopted (for comparison): worst frame |
|---|---|---|---|
| Arcsec (any binning) | 2.35″ | **2.59″** | 0.86″ |
| 1×1 pixels | 14.1 | **15.5** | 5.2 |
| 2×2 pixels | 7.1 | **7.8** | 2.6 |
| 4×4 pixels | 3.5 | **3.9** | 1.3 |

- **Where the maximum falls:** at the frame corners, 243″ from the center.
- **Its cause:** the measured rotation of 0.35–0.60° from 180°. A 0.6°
  rotation displaces the corners by 2.6″.
- **Averaged over the frame and all frames:** 1.2″.

**Assessment.** The orthogonal WCS is adequate as a first guess.

- **Compared with the pointing error:** its worst-case error (2.6″) is about
  a tenth of the error in `CRVAL` from the pointing (median 22″, up to 30″).
  The total error of the default WCS is therefore at most about 33″ in
  either case.
- **For astrometric solvers:** a 0.6° rotation is far inside the search
  ranges of astrometry.net and of `pipnick`'s matching.
- **Over time:** the rotation varies by 0.2° between nights anyway, so its
  precise value would not be stable.
- **Simplicity:** the keyword sets are simpler and easier to check, and they
  need no update if the camera's rotation changes slightly.

The adopted (rotated) matrix is about three times better at the corners. That
is worth having only if users rely on the header WCS for something more than
a starting point, e.g., to identify a target in the field before solving.

## Caveats and follow-up

- **1×1 and 4×4 are untested on sky.** They assume the plate scale is exactly
  proportional to the binning and the rotation is independent of it. Both
  hold for an ideal CCD readout. They should be checked with a few on-sky
  frames at each binning.
- **The rotation varies by about 0.2° between nights**, and should be
  monitored. If it varies more on other nights (e.g., after the camera has
  been removed and remounted), the default could be updated per run. It
  matters little for a default WCS, given the 20–30″ pointing offsets.
- **Frame of `RA`/`DEC`:**
  - The measured offsets are tens of arcseconds, not the arcminutes that
    coordinates of date would give after 26 years of precession. So `RA`/`DEC`
    are effectively J2000, and `CRVAL = RA/DEC` with `EQUINOX = 2000` is
    consistent.
  - The fit is relative to Gaia (ICRS). Writing `RADESYS = 'ICRS'` (with no
    `EQUINOX`) would be more accurate in principle, but the difference from
    FK5 J2000 (< 0.1″) is far below the pointing error.
- **Windowed frames** get no valid WCS under assumption 1. The header writer
  should omit the sky WCS for a windowed readout, or flag it, rather than
  write a wrong one.
