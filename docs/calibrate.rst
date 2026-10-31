.. _calibrate:

Photometric Calibration
=======================

If the photometry has been measured and Landolt standards were observed during
the night, pipnick can determine the photometric calibration coefficients
necessary to transform instrumental magnitudes into apparent magnitudes.  The
``pipnick_calibrate`` script finds the Landolt standard stars in the photometry
source tables and fits the photometric calibration of each band of each night.
The calibration is saved to an editable file, assessed with the residuals of the
standards, and, optionally, applied to the source tables and recorded in the
image headers.

The standards are identified by their celestial coordinates, so
``pipnick_astrometry`` (see :ref:`astrometry`) and ``pipnick_photometry`` (see
:ref:`photometry`) must be run first.

The algorithm
-------------

1. **Find the standards.**  The source tables of all of the science images are
   matched to the Landolt catalog (within ``--match_radius``).  Any image with a
   standard in its field is used for calibration.  Images whose source tables
   have no celestial coordinates are skipped with a warning.  Source tables that
   are out of date (either because their image or astrometric solution changed
   after they were measured) raise an error, reporting that the photometry needs
   to be reanalyzed using ``pipnick_photometry``.

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
       using ``--mag aper``, if this is a concern.

3. **Fit each band.**  The standards observed in each band (U, B, V, R, or I)
   during one night (UT date) are fit with the model

   .. math::

       m - M = \sum_t k_t\ x_t,

   where :math:`m = -2.5\log(f / t)` is the instrumental magnitude (for the flux
   :math:`f` in electrons and the exposure time :math:`t` in seconds), :math:`M`
   is the catalog magnitude, and the terms are chosen with ``--terms`` from the
   zero point (``zp``), the extinction (``airmass``, :math:`x = X`), the color
   term (``color``, :math:`x = C`), the quadratic color term (``color2``,
   :math:`x = C^2`), and the color dependence of the extinction
   (``airmass_color``, :math:`x = X C`).  The default is ``zp airmass color``.
   :math:`C` is the catalog color of the star; see ``--color``.  Each
   measurement is weighted by its instrumental and catalog errors and a
   systematic error (``--sys_err``), all combined in quadrature.  Outliers are
   rejected iteratively (``--rej``, ``--maxiter``).

4. **Check the fit.**  The calibration of a band fails (and the other bands
   are still calibrated) if:

   - there are fewer than ``--min_stars`` stars included in the fit,

   - the color span is less than ``--min_color_span`` while a color term is
     fit,

   - the fit is undetermined because there are as many or more terms to fit as
     there are measurements,

   - or the fitted terms are degenerate (e.g., the color does not vary).

   The ``--force`` option downgrades the first two to warnings.  Large
   coefficient errors and a negative extinction are reported as warnings.  The
   run stops with an error only if no band can be calibrated.  A night whose
   calibration fails for another reason (e.g., an unreadable source table) is
   recorded as failed in the summary and the verdict, and the other nights are
   still calibrated.

5. **Assess the fit.**  The residuals (calibrated minus catalog magnitudes)
   of every standard measurement are recorded, along with *leave-one-out*
   residuals.  For the latter, the calibrated magnitude of each star is
   predicted by refitting the calibration after excluding that star from the
   fit.  This provides an estimate of the error in the calibration that is
   particularly useful when the calibration is only based on a few stars.  The
   drift of the residuals with time is measured and reported (see below).

.. _calibrate-extinction:

Fixed terms and the nominal extinction
--------------------------------------

Any term can be fixed using the ``--fix`` option.  For example, ``--fix
color=0.03`` will fix the color term for all bands; ``--fix
airmass=B:0.28,V:0.17`` will fix the airmass terms for B and V separately.
Setting ``--fix airmass=nominal`` fixes the extinction at its nominal value for
the observatory.

.. warning::

    **If the standards of a band do not span at least** ``--min_airmass_span``
    **in airmass (0.2 by default), the extinction cannot be measured, and it is
    fixed at its nominal value.**  The ``airmass_color`` term, if used, is also
    fixed at 0.  This *extinction fallback* is expected to happen often: on most
    nights, the standards are observed over a small range of airmass.

    The "nominal" extinction is not a measurement for the night.  It comes from
    an extinction curve for Mt. Hamilton integrated over generic filter response
    curves (see below); this means the calibration may have systematic errors
    that are not included in the error budget.

    When this fallback is invoked, it is recorded in every product, so that it
    cannot be missed:

    - in the calibration file, the ``source`` of the ``airmass`` term is set to
      ``"fallback"``, with a ``fallback_reason``, under a comment stating that
      the extinction was **not measured**,

    - ``PNCALEXT = 'fallback'`` is added to the image headers,

    - the ``extinction`` metadata is included in the residual table,

    - an "EXTINCTION NOT MEASURED" label is added to every panel of the band in
      the relevant QA figures,

    - the ``extinction`` column of the summary table is ``fallback``,

    - a banner at the end of the log lists every band that fell back.

    To avoid the fallback, observe standards over a wider range of airmass.  To
    make the band fail instead, use ``--no_fallback``.

The nominal extinction of each band is computed by
:mod:`pipnick.photometry.extinction` from the extinction curve of the
observatory, named by the ``extinction`` entry of the ``[observatory]`` section
of the camera configuration.  For Lick Observatory, this is
``pipnick/data/extinction/mthamextinct.dat``, copied from `PypeIt
<https://github.com/pypeit/PypeIt>`__, which covers 3200 to 10870 angstroms.
Its origin is not documented (the PypeIt list of extinction curves notes that
they were "Generally grabbed from the Internet.  Buyer beware.").  The curve is
integrated over the filter response curves assuming a source with a flat
spectrum in :math:`f_\lambda`, giving (in magnitudes per airmass) U 0.604, B
0.300, V 0.178, R 0.120, and I 0.075.  The curve covers only part of the U band,
so U is the least certain.

.. note::

    The filter response curves are the generic Bessell UBVRI curves provided by
    the ``speclite`` package, not those specific to the Nickel filters.  They
    can be replaced by Lick/Nickel-specific curves when these are available.

Uncertainties
-------------

All of the reported errors --- the coefficients, the zero points, and the
calibrated magnitudes --- are **statistical only**.  Fixed and fallback terms
have no error, and the systematic error they introduce is not included.  The
leave-one-out residuals can provide a better estimate of the actual accuracy of a
calibration from few stars.

Drift of the zero point
-----------------------

The calibration assumes that the zero point is constant through the night.  The
drift of the residuals with time is measured (per band, and for all bands
together) and recorded in the log, the calibration file, the summary, and the QA
figures.  A warning is issued if it is significant and larger than
``--drift_warn`` (mag/hour).  The drift is *not* fit or included in the
photometric calibration determination.

Inputs
------

- **The ECSV metadata table** (``raw_table``), written by ``pipnick_metadata``;
  see :ref:`metadata`.  It provides the airmass and observation time of each
  image.

  One run calibrates the images of one table, i.e., one instrument setup.  If
  the table includes frames from more than one night (UT date), each night is
  fit separately; no fit combines nights.  Nights and setups in other tables
  are calibrated by running ``pipnick_calibrate`` on each table.

- **The photometry source tables** of the science images, with celestial
  coordinates (see :ref:`photometry` and :ref:`astrometry`).  They must be up to
  date.

- **Images of Landolt standard fields** taken during the night, in U, B, V, R,
  and/or I.  The more standards, and the wider their range of airmass and color,
  the better the calibration.

- **The Landolt catalog**, distributed with pipnick (see below).

Execution
---------

An example call looks like this:

.. code-block:: console

   pipnick_calibrate n1_2x2_slow_full.ecsv --apply both

The required argument is the per-setup metadata table.

The options are:

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
  0.05).  Standards without a catalog error are used, but a warning is issued.

- ``--min_cat_obs``: The minimum number of Landolt observations of a standard
  (default 1).

- ``--sys_err``: The systematic error (mag) added to each measurement (default
  0.01).

- ``--rej``, ``--maxiter``: The rejection threshold in standard deviations
  (default 3) and the maximum number of rejection iterations (default 5) for
  excluding observations from the calibration fit.

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

Observations taken using filters other than UBVRI (e.g., SDSS filters) currently
cannot be calibrated with the Landolt standards; their images are skipped with a
warning.

Outputs
-------

The output files are added to the ``calibrate`` subdirectory of the reduction
directory:

- ``<setup>_<night>.toml`` provides the calibration results (see below).

- ``<setup>_<night>_residuals.ecsv`` is a table with one row per matched
  standard measurement.  Each row provides the image, band, star, match
  separation, airmass, MJD, catalog magnitude and color, instrumental magnitude,
  the number of residual detections near the standard (``n_residual``), the
  calibrated magnitude (``mag_cal``), the residual (``resid``, calibrated minus
  catalog) and its error, the leave-one-out residual (``resid_loo``; masked if
  it cannot be computed), whether it was used or rejected, and the reason it is
  unusable.

- ``summary.ecsv`` provides the calibration status of each band, the numbers of
  stars and measurements, the rms and leave-one-out rms of the residuals, how
  the extinction was set (``fit``, ``fixed``, ``fallback``, ``none`` if there is
  no airmass term, or ``failed``), the extinction coefficient, the zero point
  above the atmosphere (airmass 0, for color 0) and its error (``zp``,
  ``zp_err``), the drift, and what was done with the calibration file
  (``action``: ``fit``, ``up to date``, ``edited (kept)``, or ``failed``).

- ``qa/<setup>_<night>_residuals.png`` shows the residuals of each band against
  airmass, color, catalog magnitude, and time.

- ``<image>_cal.ecsv`` provides the calibrated source table
  of each science image, if ``--apply catalog`` is set.

If ``--apply image`` is set, the calibration keywords are also added to the
headers of the reduced images; see :ref:`calibrate-applied`.

The log is ``log/pipnick_calibrate_<time>.log`` and placed in the reduction
directory; see :ref:`common`.

The calibration file
--------------------

The calibration is saved as a TOML file so that it can be read and edited.  Its
``[calibration]`` table records the setup, night, camera, catalog, the
instrumental magnitude calibrated, the creation time, the pipnick version
(``pipnick_version``), the overall drift, and three hashes.  The
``catalog_hash`` is for the standard-star catalog, the ``input_hash`` is for the
source tables and settings used, and the ``fit_hash`` is for the content of the
file.  Each band has a ``[band.<band>]``
table with its ``status`` (``ok`` or ``failed`` with a ``reason``), the color,
the ``terms``, their ``coefficients``, ``errors``, ``source`` (``fit``,
``fixed``, or ``fallback``), and ``covariance``.  Each band also provides an
assessment of the results providing the numbers of stars, measurements, rejected
and unusable measurements, the rms, reduced chi-square, and leave-one-out rms of
the residuals, the drift, the airmass and color ranges, and the images used.

All of this may be edited.  Note, however, when ``pipnick_calibrate``
is rerun:

- if the inputs are unchanged, the calibration is not refit;

- if the inputs changed, the night is refit;

- if the file has been edited (it no longer matches its ``fit_hash``), it is
  kept, and used by ``--apply``, with a warning if its inputs have changed;

- with ``--overwrite``, the night is always refit.

.. _calibrate-applied:

Applying the calibration
------------------------

The calibration determined from the photometric standards can be used to
calibrate the instrumental magnitudes for other sources obtained during the
night.  When setting ``--apply catalog``, the source tables of all science
images are calibrated.  The sources in the images of each object are identified
across the images by their celestial coordinates, and the calibration is
iterated to find the color of each source from its calibrated magnitudes,
averaged over the images in each band.  A source without a measured color (e.g.,
observed in only one band) uses the default color (``--default_color``) and is
flagged.

The calibrated source tables, ``calibrate/<image>_cal.ecsv``, have the source
``id``, ``x``, ``y``, ``ra``, ``dec``, the instrumental magnitude, the
calibrated magnitude (``mag``) and its statistical error (``mag_err``),
``color_used``, ``color_assumed``, and the identifier of the source across the
images (``source_id``).  The original source tables with instrumental magnitudes
from the photometry itself are not changed.

With ``--apply image``, the calibration is recorded in the headers of the
reduced images; their pixels are not changed.  The calibrated magnitude of a
source with count rate :math:`r` (electrons per second) and color :math:`C` is

.. math::

    M = -2.5\log(r) + {\rm PNCALZP} + {\rm PNCALCT}\ C + {\rm PNCALC2}\ C^2,

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

- ``PNCALVER``: the pipnick version that applied the calibration;

- ``PHOTZP``: the zero point for :math:`C = 0`, equal to ``PNCALZP``.

The keywords are rewritten only when the calibration file changes.  The keywords
are removed from images whose band is not calibrated.  If an image is reduced
again, its keywords are lost, meaning you will need to rerun ``pipnick_calibrate
--apply image`` to restore them from the calibration file.

The Landolt catalog
-------------------

The catalog distributed with pipnick (``pipnick/data/catalogs/landolt.ecsv``)
merges Landolt (2009, AJ, 137, 4186; VizieR ``J/AJ/137/4186``) and Landolt
(2013, AJ, 146, 131; VizieR ``J/AJ/146/131``) and provides 942 stars with UBVRI
photometry, accurate coordinates, and proper motions.  The positions are
propagated to the epoch of each observation.  The UBRI magnitudes are derived
from the published V magnitudes and colors.  The original Landolt (1992) catalog
(VizieR ``II/183A``) is not used: its coordinates are only accurate to about 4
arcsec and it has no proper motions.  The catalog was constructed using
:func:`~pipnick.photometry.landolt.build_catalog`.

Note that some standards (e.g., in the PG1530+057 field) were observed only once
by Landolt and have no published errors.  They are used by default; however, you
can exclude them by setting ``--min_cat_obs 2``.

This research has made use of the VizieR catalogue access tool, CDS,
Strasbourg, France (Ochsenbein et al. 2000, A&AS, 143, 23).

Troubleshooting
---------------

- **"No source table has celestial coordinates"**:  The standards are found by
  their coordinates, which means that ``pipnick_astrometry`` and
  ``pipnick_photometry`` must have been executed (in that order) to ensure
  coordinates are included in the source tables.

- **A source table "is out of date"**: Rerun ``pipnick_photometry`` to update
  it.

- **A band fails with too few stars**: A night with few standards cannot
  constrain every term.  To obtain a solution, try fitting fewer terms (e.g.,
  set ``--terms zp`` with the extinction fixed at its nominal value by setting
  ``--fix airmass=nominal``), or lower the ``--min_stars`` value and accept a
  less certain calibration.  You should check the leave-one-out residuals to
  determine if the calibration is viable.

- **A band fails on its color span**: The standards have similar colors, meaning
  that a color term is too uncertain to be included in the fit.  Fix it (``--fix
  color=<value>``) or remove it (``--terms``).

- **The extinction was not measured**: See :ref:`calibrate-extinction`.

- **"No band could be calibrated"**:  The run stops with an error due to reasons
  provided in the summary and the calibration file.

- **A band that cannot be calibrated**:  This may be because no standards were
  observed in it.  The code will not calibrate the relevant images, and
  ``pipnick_status`` will continue to suggest running ``pipnick_calibrate`` for
  them.  This is currently an unsolved issue; see :ref:`issues`.

- **Filters other than UBVRI**: Unfortunately, these cannot yet be calibrated
  with pipnick.  We will support use of other standard catalogs for photometric
  calibration in the future.

- **The PSF magnitudes of bright standards** are less precise than their
  aperture magnitudes; consider using ``--mag aper`` (see the note above).

Using the pipeline in Python
----------------------------

.. code-block:: python

   from pipnick.pipelines.calibrate import calibrate_all

   files = calibrate_all('n1_2x2_slow_full.ecsv', apply='both')

See :func:`~pipnick.pipelines.calibrate.calibrate_all` for the full set of
parameters, and :mod:`pipnick.photometry.calibrate` for the functions used to
fit and apply the calibration.
