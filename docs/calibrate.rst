.. _calibrate:

pipnick_calibrate: photometric calibration
==========================================

``pipnick_calibrate`` finds the Landolt standard stars in the photometry
source tables and fits the photometric calibration of each band of each
night.  The calibration is saved to an editable file, assessed with the
residuals of the standards, and, optionally, applied to the source tables and
recorded in the image headers.

The standards are identified by their celestial coordinates, so
``pipnick_astrometry`` (see :doc:`astrometry`) and ``pipnick_photometry`` (see
:doc:`photometry`) must be run first.

The algorithm
-------------

1. **Find the standards.**  The source tables of all of the science images
   are matched to the Landolt catalog (within ``--match_radius``), so any
   image with a standard in its field is used, whatever its object.  Images
   whose source tables have no celestial coordinates are skipped with a
   warning.  A source table that is out of date (its image or astrometric
   solution changed after it was measured) is an error: rerun
   ``pipnick_photometry``.

2. **Select the measurements.**  A matched standard is not used if its
   catalog magnitude or color is missing, its instrumental magnitude is
   masked, its fit is flagged by ``photutils``, its aperture includes masked
   pixels (e.g., saturated pixels or a bad column), it is contaminated by a
   neighbor (see below), its catalog error is larger than ``--max_cat_err``,
   or it was observed fewer than ``--min_cat_obs`` times by Landolt.  The
   reason is recorded for every measurement.

   The PSF photometry of bright stars often detects faint, spurious sources in
   the residuals of their fits, within about one FWHM of the star, which are
   then fit together with it.  These *residual detections* (sources detected
   in a later iteration of the PSF photometry, within ``--residual_radius``
   FWHM) do not exclude a standard; their number is recorded
   (``n_residual``).  A standard is *contaminated* if any other source within
   its photometric aperture has more than ``--max_neighbor_flux`` (1% by
   default) of its flux.

   .. note::

       The residual detections take a few percent of the flux of the star in
       the PSF fit, so the PSF magnitudes of bright standards are less precise
       than their aperture magnitudes.  In back-to-back frames of the
       PG1530+057 standards (2024-06-24), the PSF magnitudes differed by up to
       0.045 mag, and the aperture magnitudes by at most 0.012 mag.  Consider
       ``--mag aper`` (and aperture photometry of the science sources).

3. **Fit each band.**  The standards observed in each band (U, B, V, R, or I)
   during one night (UT date) are fit with the model

   .. math::

       m - M = \sum_t k_t\ x_t,

   where :math:`m = -2.5\log_{10}(f / t)` is the instrumental magnitude (for
   the flux :math:`f` in electrons and the exposure time :math:`t`),
   :math:`M` is the catalog magnitude, and the terms are chosen with
   ``--terms`` from the zero point (``zp``), the extinction (``airmass``,
   :math:`x = X`), the color term (``color``, :math:`x = C`), the quadratic
   color term (``color2``, :math:`x = C^2`), and the color dependence of the
   extinction (``airmass_color``, :math:`x = X C`).  The default is ``zp
   airmass color``.  :math:`C` is the catalog color of the star; see
   ``--color``.  Each measurement is weighted by its instrumental and catalog
   errors and a systematic error (``--sys_err``), added in quadrature.
   Outliers are rejected iteratively (``--rej``, ``--maxiter``).

4. **Check the fit.**  The calibration of a band fails (and the other bands
   are still calibrated) if:

   - there are fewer than ``--min_stars`` stars;
   - the color span is less than ``--min_color_span`` while a color term is
     fit;
   - there are no more measurements than fitted terms, or the fitted terms
     are degenerate (e.g., the color does not vary).

   ``--force`` downgrades the first two to warnings.  Large coefficient
   errors and a negative extinction are reported as warnings.  The run stops
   with an error only if no band can be calibrated.  A night whose
   calibration fails for another reason (e.g., an unreadable source table) is
   recorded as failed in the summary and the verdict, and the other nights
   are still calibrated.

5. **Assess the fit.**  The residuals (calibrated minus catalog magnitudes)
   of every standard measurement are recorded, along with *leave-one-out*
   residuals: each star is predicted by the calibration fit without it, which
   is an honest estimate of the error when there are few stars.  The drift of
   the residuals with time is measured and reported (see below).

Fixed terms and the nominal extinction
--------------------------------------

Any term can be fixed with ``--fix``, e.g., ``--fix color=0.03`` (for all
bands) or ``--fix airmass=B:0.28,V:0.17``.  ``--fix airmass=nominal`` fixes the
extinction at its nominal value for the observatory.

.. warning::

    **If the standards of a band do not span at least ``--min_airmass_span``
    in airmass (0.2 by default), the extinction cannot be measured, and it is
    fixed at its nominal value instead** (and the ``airmass_color`` term, if
    used, is fixed at 0).  This *extinction fallback* is common: on most
    nights, the standards are observed over a small range of airmass.

    The nominal extinction is not a measurement for the night.  It comes from
    an extinction curve for Mt. Hamilton folded through generic filter
    responses (see below), so the calibration may have systematic errors that
    are not included in its errors.

    The fallback is recorded in every product, so that it cannot be missed:

    - the calibration file: the ``source`` of the ``airmass`` term is
      ``"fallback"``, with a ``fallback_reason``, under a comment stating that
      the extinction was **not measured**;
    - the image headers: ``PNCALEXT = 'fallback'``;
    - the residual table: its ``extinction`` metadata;
    - the QA figures: "EXTINCTION NOT MEASURED" on every panel of the band;
    - the summary: its ``extinction`` column, and a banner at the end of the
      output listing every band that fell back.

    To avoid the fallback, observe standards over a wider range of airmass
    (e.g., the same field early and late in the night).  To make the band
    fail instead, use ``--no_fallback``.

The nominal extinction of each band is computed by
:mod:`pipnick.photometry.extinction` from the extinction curve of the
observatory, named by the ``extinction`` entry of the ``[observatory]``
section of the camera configuration.  For Lick Observatory, this is
``pipnick/data/extinction/mthamextinct.dat``, copied from `PypeIt
<https://github.com/pypeit/PypeIt>`__, which covers 3200 to 10870 angstroms.
Its origin is not documented (the PypeIt list of extinction curves notes that
they were "Generally grabbed from the Internet.  Buyer beware.").  The curve
is folded through the filter response for a source with a flat spectrum in
:math:`f_\lambda`, giving (in magnitudes per airmass) U 0.604, B 0.300, V
0.178, R 0.120, and I 0.075.  The curve covers only part of the U band, so U
is the least certain.

.. note::

    The filter responses are the generic Bessell UBVRI curves provided by the
    ``speclite`` package, not those of the Nickel filters.  They can be
    replaced by Lick/Nickel-specific curves when these are available.

The errors
----------

All of the reported errors (of the coefficients, the zero points, and the
calibrated magnitudes) are **statistical only**.  Fixed and fallback terms
have no error, and the systematic error they introduce is not included.  The
leave-one-out residuals are a better guide to the actual accuracy of a
calibration from few stars.

Drift of the zero point
-----------------------

The calibration assumes that the zero point is constant through the night.
The drift of the residuals with time is measured (per band, and for all
bands together), reported in the log, the calibration file, the summary, and
the QA figures, and a warning is issued if it is significant and larger than
``--drift_warn`` (mag/hour).  The drift is *not* fit; adding a time term to the
model is an area for improvement.

Inputs
------

- **A per-setup metadata table** (``-t``), used to reduce the images (see
  :doc:`metadata` and :doc:`reduce`); it sets the airmass and observation
  time of each image.
- **The source tables** of the science images, with celestial coordinates
  (see :doc:`photometry` and :doc:`astrometry`).  They must be up to date.
- **Images of Landolt standard fields** taken during the night, in U, B, V, R,
  or I.  The more standards, and the wider their range of airmass and color,
  the better the calibration.
- **The Landolt catalog**, distributed with pipnick (see below).

Running it
----------

.. code-block:: console

   pipnick_calibrate -t n1_2x2_slow_full.ecsv --apply both

The options are (the logging options are common to every script; see
:doc:`common`):

- ``-t``, ``--raw_table``: The per-setup metadata file used to reduce the
  images (required).
- ``-r``, ``--rdxdir``: The reduction directory.  If not provided, the
  directory saved in the metadata file is used.
- ``--images``: Only search these images for standards, given by file name
  (with or without the extension) or combination group (``comb_id``).  This
  does not limit the images the calibration is applied to.
- ``--terms``: The terms of the model (default ``zp airmass color``).
- ``--fix``: Fix terms; see above.
- ``--color``: The color used to calibrate a band, as ``<band>=<color>``, e.g.
  ``R=R-I``.  The defaults are U-B for U, B-V for B and V, V-R for R, and V-I
  for I.
- ``--mag``: The instrumental magnitudes to calibrate: ``psf`` (default; the
  ``mag_inst`` column) or ``aper`` (from the aperture fluxes).  No aperture
  correction is applied, so a calibration of aperture magnitudes is valid only
  for magnitudes measured in the same aperture; a warning is issued.
- ``--match_radius``: The radius (arcsec) for matching sources to the
  standards (default 2).
- ``--max_neighbor_flux``: The fraction of the flux of a standard above which
  a neighbor within its aperture contaminates it (default 0.01).
- ``--residual_radius``: The radius (FWHM) within which sources detected in the
  residuals of a standard do not exclude it (default 1.5).
- ``--max_cat_err``: The maximum catalog error (mag) of a standard (default
  0.05).  Standards without a catalog error are used, with a warning.
- ``--min_cat_obs``: The minimum number of Landolt observations of a standard
  (default 1).
- ``--sys_err``: The systematic error (mag) added to each measurement (default
  0.01).
- ``--rej``, ``--maxiter``: The rejection threshold in standard deviations
  (default 3) and the maximum number of rejection iterations (default 5).
- ``--min_stars``, ``--min_airmass_span``, ``--min_color_span``: The sample
  checks (defaults 3, 0.2, and 0.3 mag).
- ``--force``: Downgrade the minimum-stars and color-span errors to warnings.
- ``--no_fallback``: Fail instead of using the nominal extinction.
- ``--drift_warn``: The drift (mag/hour) above which a significant drift is
  reported as a warning (default 0.02).
- ``--no_qa``: Do not write the QA figures to ``calibrate/qa``.  By default, a
  figure is written when a night is fit, and drawn from the residual table for
  a night that is not refit but lacks one.
- ``-o``, ``--overwrite``: Refit, even if the inputs have not changed or the
  calibration file has been edited.
- ``--apply``: Apply the calibration to the science images of the night:
  ``catalog``, ``image``, or ``both``; see below.
- ``--default_color``: The color assumed for sources without a measured color,
  as ``<color>=<value>``, e.g. ``B-V=0.6`` (default 0).

Filters other than U, B, V, R, and I (e.g., the primed filters of the Andor
camera) cannot be calibrated with the Landolt standards; their images are
skipped with a warning.

Outputs
-------

The products are in the ``calibrate`` subdirectory of the reduction directory:

- ``<setup>_<night>.toml``: the calibration (see below);
- ``<setup>_<night>_residuals.ecsv``: one row per matched standard
  measurement, with the image, band, star, match separation, airmass, MJD,
  catalog magnitude and color, instrumental magnitude, the number of residual
  detections near the standard (``n_residual``), the calibrated
  magnitude (``mag_cal``), the residual (``resid``, calibrated minus catalog)
  and its error, the leave-one-out residual (``resid_loo``; masked if it
  cannot be computed), whether it was used or rejected, and the reason it is
  unusable (``-`` if it is usable);
- ``summary.ecsv``, also logged: the status of each band, the numbers of stars
  and measurements, the rms and leave-one-out rms of the residuals, how the
  extinction was set (``fit``, ``fixed``, ``fallback``, or ``failed``), the
  extinction coefficient, the zero point above the atmosphere (airmass 0, for
  color 0) and its error (``zp``, ``zp_err``), the drift, and what was done
  with the calibration file (``action``: ``fit``, ``up to date``, or
  ``edited (kept)``);
- ``qa/<setup>_<night>_residuals.png`` (unless ``--no_qa``): the residuals of
  each band against airmass, color, catalog magnitude, and time;
- with ``--apply catalog``, ``<image>_cal.ecsv``: the calibrated source table
  of each science image; and with ``--apply image``, the calibration keywords
  in the headers of the reduced images (see `Applying the calibration`_).

The log is ``log/pipnick_calibrate_<time>.log`` in the reduction directory;
see :doc:`common`.

The calibration file
--------------------

The calibration is saved as TOML, so that it can be read and edited.  Its
``[calibration]`` table records the setup, night, camera, catalog, the
instrumental magnitude calibrated, the creation time, the overall drift, and
two hashes: ``input_hash`` (of the source tables and settings used) and
``fit_hash`` (of the content of the file).  Each band has a ``[band.<band>]``
table with its ``status`` (``ok`` or ``failed``, with a ``reason``), the
color, the ``terms``, their ``coefficients``, ``errors``, ``source``
(``fit``, ``fixed``, or ``fallback``), and ``covariance``, and the
assessment: the numbers of stars, measurements, rejected and unusable
measurements, the rms, reduced chi-square, and leave-one-out rms of the
residuals, the drift, the airmass and color ranges, and the images used.

The coefficients (and the terms) may be edited.  When ``pipnick_calibrate``
is rerun:

- if the inputs are unchanged, the calibration is not refit;
- if the inputs changed, the night is refit;
- if the file has been edited (it no longer matches its ``fit_hash``), it is
  kept, and used by ``--apply``, with a warning if its inputs have changed;
- with ``--overwrite``, the night is always refit.

Applying the calibration
------------------------

With ``--apply catalog``, the source tables of the science images are
calibrated.  The sources in the images of each object are identified across
the images by their celestial coordinates, and the calibration is iterated to
find the color of each source from its calibrated magnitudes, averaged over
the images in each band.  A source without a measured color (e.g., observed
in only one band) uses the default color (``--default_color``) and is flagged.
The calibrated tables, ``calibrate/<image>_cal.ecsv``, have the source
``id``, ``x``, ``y``, ``ra``, ``dec``, the instrumental magnitude, the
calibrated magnitude (``mag``) and its statistical error (``mag_err``),
``color_used``, ``color_assumed``, and the identifier of the source across
the images (``source_id``).  The source tables of the photometry are not
changed.

With ``--apply image``, the calibration is recorded in the headers of the
reduced images; their pixels are not changed.  The calibrated magnitude of a
source with count rate :math:`r` (electrons per second) and color :math:`C` is

.. math::

    M = -2.5\log_{10}(r) + {\rm PNCALZP} + {\rm PNCALCT}\ C
        + {\rm PNCALC2}\ C^2,

with the keywords:

- ``PNCALZP``, ``PNCALZE``: the zero point at the airmass of the image, and
  its statistical error;
- ``PNCALCT``, ``PNCALC2``: the color coefficients;
- ``PNCALCOL``: the color, e.g. ``B-V``;
- ``PNCALMAG``: the instrumental magnitude calibrated (``psf`` or ``aper``);
- ``PNCALSYS``: the magnitude system (``Vega``);
- ``PNCALEXT``: how the extinction was set: ``fit``, ``fixed``, ``fallback``
  (the nominal value; not measured), or ``none``;
- ``PNCALFIL``, ``PNCALHSH``: the calibration file and its hash;
- ``PHOTZP``: the zero point for :math:`C = 0`, equal to ``PNCALZP``.

The keywords are rewritten only when the calibration file changes, and are
removed from images whose band is not calibrated.  They do not cause the
photometry to be measured again.  If an image is reduced again, its keywords
are lost; rerun ``pipnick_calibrate --apply image`` to restore them from the
calibration file.

The Landolt catalog
-------------------

The catalog distributed with ``pipnick``, ``pipnick/data/catalogs/landolt.ecsv``,
merges Landolt (2009, AJ, 137, 4186; VizieR ``J/AJ/137/4186``) and Landolt
(2013, AJ, 146, 131; VizieR ``J/AJ/146/131``): 942 stars with UBVRI
photometry, accurate coordinates, and proper motions.  The positions are
propagated to the epoch of each observation.  The U, B, R, and I magnitudes
are derived from the published V magnitudes and colors.  The original Landolt
(1992) catalog (VizieR ``II/183A``) is not used: its coordinates are only
accurate to about 4 arcsec, and it has no proper motions.  The catalog was
retrieved with :func:`~pipnick.photometry.landolt.build_catalog`.

Note that some standards (e.g., in the PG1530+057 field) were observed only
once by Landolt, and have no published errors.  They are used by default; use
``--min_cat_obs 2`` to exclude them.

This research has made use of the VizieR catalogue access tool, CDS,
Strasbourg, France (Ochsenbein et al. 2000, A&AS, 143, 23).

Troubleshooting
---------------

- **"No source table has celestial coordinates"**: the standards are found by
  their coordinates.  Run ``pipnick_astrometry``, then ``pipnick_photometry``
  (which adds the coordinates).
- **A source table "is out of date"**: rerun ``pipnick_photometry``.
- **A band fails with too few stars** (``--min_stars``, 3 by default): a
  night with few standards cannot constrain every term.  Fit fewer terms
  (e.g., ``--terms zp``, with the extinction fixed at its nominal value by
  ``--fix airmass=nominal`` if it is in the model), or lower ``--min_stars``
  and accept a less certain calibration (check its leave-one-out rms).
- **A band fails on its color span**: the standards have similar colors, so a
  color term cannot be fit.  Fix it (``--fix color=<value>``) or remove it
  (``--terms``).
- **The extinction was not measured** (the banner at the end of the log):
  see the warning above.
- **"No band could be calibrated"**: the run stops with an error; the reasons
  are in the summary and the calibration file.
- **A band that cannot be calibrated** (e.g., no standards observed in it):
  its images are not calibrated, and ``pipnick_status`` keeps suggesting
  ``pipnick_calibrate`` for them; see :doc:`issues`.
- **Filters other than UBVRI** (e.g., the primed filters of the Andor camera)
  cannot be calibrated with the Landolt standards.
- **The PSF magnitudes of bright standards** are less precise than their
  aperture magnitudes; consider ``--mag aper`` (see the note above).

Using the pipeline in Python
----------------------------

.. code-block:: python

   from pipnick.pipelines.calibrate import calibrate_all

   files = calibrate_all('n1_2x2_slow_full.ecsv', apply='both')

See :func:`~pipnick.pipelines.calibrate.calibrate_all` for the full set of
parameters, and :mod:`pipnick.photometry.calibrate` for the functions used to
fit and apply the calibration.
