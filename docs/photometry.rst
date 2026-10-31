.. _photometry:

Point-Spread-Function (PSF) and Aperture Photometry
===================================================

After performing basic image processing, pipnick can perform point-source
photometry using ``pipnick_photometry``.  This script detects and measures the
flux for sources in each reduced science image with point-spread-function (PSF)
and aperture photometry.  The results are written to source tables for each
image.  It does not require an astrometric solution: the source coordinates can
be provided in pixels only.  However, given that most of the subsequent
processing steps require the astrometric solution, we encourage users to perform
the astrometry (see :doc:`astrometry`) before the photometry.

.. _photometry-algorithm:

The algorithm
-------------

For each image, photometry proceeds as follows:

1. **Model the PSF.**  To perform PSF photometry, the code must first build an
   empirical measurement of the point-spread function using isolated bright
   stars.  These stars are detected, their local background is subtracted, and
   they are stacked to form a high S/N PSF measurement, which is then modeled by
   a circular or elliptical Moffat function.  The fit sets the FWHM used for all
   the scales below.

2. **Subtract the background.**  A smooth, 2D background is subtracted
   from the full image.

3. **PSF photometry.**  Sources are detected using
   :class:`photutils.detection.DAOStarFinder`, identifying all sources above
   ``--thresh`` times the background noise.  Sources that are no closer than
   ``--min_sep`` FWHM are independently modeled by the (fixed-shape) PSF model,
   whereas sources closer than ``--group_sep`` FWHM are fit simultaneously.  The
   fit is performed iteratively (using
   :class:`photutils.psf.IterativePSFPhotometry`; see ``--mode``) where each
   iteration fits and subtracts the model sources, and the subsequent iteration
   detects new sources in the residual image.  With ``--group``, sources fit
   together are consolidated into a single source.

4. **Vet the sources.**  Sources are removed if their fit failed, their fitted
   positions that are either outside the image or moved by more than one FWHM
   from where their original detected location, or their fitted fluxes are not
   positive.  These soucres are typically spurious detections in the residuals
   of bright stars.  They are listed, with the reason, in
   ``photometry/<image>_removed.ecsv``.

5. **Aperture photometry** is also performed by summing the
   background-subtracted flux in a circular aperture of radius ``--aper_radius``
   times the FWHM.

The image uncertainties used to weight the fits and compute the flux errors are
the uncertainties propagated by ``pipnick_reduce`` (the ``UNCERT`` extension of
the reduced image).  This includes the read noise, the Poisson noise of the
counts,  and all uncertainties propagated through the image processing steps
(overscan and bias subtraction, and field-flattening).

Inputs
------

- **The ECSV metadata table** (``raw_table``), written by ``pipnick_metadata``;
  see :doc:`metadata`.  The pointing of each image (``ra``, ``dec``) is taken
  from it.

- **The reduced science images**, with their uncertainties (``UNCERT``) and
  masks (``MASK``).

- Optionally, **the astrometric solution** (see :doc:`astrometry`), for the
  celestial coordinates of the sources.  The photometry can be run before or
  after the astrometry: rerunning the photometry after obtaining an astrometric
  solution only adds the coordinates to the output source catalog.

Execution
---------

An example call looks like this:

.. code-block:: console

   pipnick_photometry n1_2x2_slow_full.ecsv

The required argument is the per-setup metadata table.

The options are:

- ``-r``, ``--rdxdir``: The reduction directory.  If not provided, the
  directory saved in the metadata file is used.

- ``--images``: A restricted set of images to analyze, given by file name (with
  or without the extension) or combination group (``comb_id``).

- ``-o``, ``--overwrite``: Reanalyze the images, even if they are up to date.

- ``--no_qa``: Do not write the QA figures to ``photometry/qa``.  By default,
  the figures are drawn for each measured image, including for images that were
  previously processed but have no QA figures.

- ``--thresh``: The source detection threshold, as a multiple of the
  background noise (default 8).

- ``--psf_thresh``: The detection threshold for the stars used to model the
  PSF (default 15).

- ``--fwhm``: The initial guess for the FWHM in arcsec (default: typical for
  the camera).

- ``--fittype``: The PSF model shape to use.  Must be ``circ`` (default) or
  ``ellip``.

- ``--stack``: The method used to stack the PSF stars.  Must be ``mean``
  (default) or ``median``.

- ``--mode``: The iterative fitting mode.  Must be ``all`` (default) or ``new``;
  see :class:`photutils.psf.IterativePSFPhotometry`.

- ``--min_sep``, ``--group_sep``: The minimum separation of detected sources,
  and the separation below which sources are fit together.  Both are in units of
  the FWHM.

- ``--group``: Consolidate the sources fit together into single sources.

- ``--aper_radius``: The radius of the photometric apertures in FWHM.

The execution time roughly scales with the number of sources identified.
However, rerunning the script is cheap: an image is measured again only if it
changed (e.g., it was reduced again) or the settings changed.  If only its
astrometric solution changed, for example because ``pipnick_astrometry`` was run
after ``pipnick_photometry``, only the celestial coordinates are updated.

Outputs
-------

A source catalog is provided for each image with the name
``photometry/<image>_phot.ecsv`` and saved to the reduction directory.  The
source catalog contains the following columns:

- ``id``, ``group_id``, ``group_size``, ``iter_detected``: The source and group
  identifiers, the number of sources fit together, and the iteration in which
  each source was detected.

- ``x``, ``y``, ``x_err``, ``y_err``: The 0-based pixel coordinates.

- ``ra``, ``dec``: The ICRS coordinates in degrees.  These are only provided if
  the image has an astrometric solution from ``pipnick_astrometry`` and they are
  masked otherwise.  The nominal WCS that the reduced images inherit from the
  raw headers is never used.

- ``flux_psf``, ``flux_psf_err``: The PSF flux in electrons.

- ``flux_aper``, ``flux_aper_err``: The aperture flux in electrons, masked if
  the aperture includes masked pixels or extends beyond the image.

- ``mag_inst``, ``mag_inst_err``: The instrumental magnitude,
  :math:`-2.5 \log(f / t)` for the PSF flux :math:`f` (``flux_psf``) and the
  exposure time :math:`t`.

- ``local_bkg``, ``qfit``, ``cfit``, ``reduced_chi2``, ``flags``: The local
  background and the fit-quality measures and flags provided by ``photutils``.

The table metadata records the image name, its object, filter, exposure time,
airmass, and MJD.  It also records the PSF model parameters, (type, parameters,
FWHM in pixels and arcsec, and the number of stars used) the number of sources
removed by the vetting, the settings used to perform the photometry, and the
provenance of the image and its astrometric solution.

All sources removed (see :ref:`photometry-algorithm`, step 4) are written to
``photometry/<image>_removed.ecsv``, with the 0-based pixel coordinates where
they were detected (``x_init``, ``y_init``), their fitted coordinates and flux
(``x_fit``, ``y_fit``, ``flux_fit``, which are masked if the fit failed), the
iteration of the PSF photometry that detected them (``iter_detected``), and the
reason they were removed (``failed``, ``outside``, ``moved``, or
``nonpositive``).

A summary of the photometry is provided by ``photometry/summary.ecsv``, which
summarizes the results for each image.  It provides the photometry status
(``measured``, ``up to date``, ``coordinates updated``, ``missing``, or
``failed``, with the ``reason``), the number of sources, the FWHM in pixels and
arcsec, the number of PSF stars, and the astrometric solution method (``none``
if the coordinates are masked).  The FWHM of each image is checked against the
camera's expected range; see :doc:`common`.

There are two QA figures for each image:

- ``photometry/qa/<image>_psf.png`` shows the stacked PSF stars, the model, the
residuals, and the radial profile

- ``photometry/qa/<image>_sources.png`` shows the image with the measured
sources numbered by decreasing PSF flux.

The log is written to ``log/pipnick_photometry_<time>.log``; see :doc:`common`.

The instrumental magnitudes are *not* calibrated.  Use ``pipnick_calibrate`` to
perform photometric calibration; currently this can only be done if Landolt
standard stars were observed during the night (see :doc:`calibrate`).

Troubleshooting
---------------

- **"No sources are good enough to model the PSF"**.  This occurs if their are
  too few bright, isolated, unsaturated stars.  To successfully analyze this
  image, try lowering ``--psf_thresh`` or giving a better guess for the
  ``--fwhm``.

- **The FWHM is outside its expected range**.  This may occur if the focus or
  tracking was bad.  Inspect the ``*_psf.png`` figure and potentially skip this
  image or try fitting an elongated PSF (``--fittype ellip``).

- **No celestial coordinates** (``wcs`` is ``none`` in the summary): the image
  has no astrometric solution.  Run ``pipnick_astrometry``, then rerun the
  photometry, which only adds the coordinates.

- **Many removed sources around bright stars**.  The PSF fits leave residuals
  around bright stars, which are then detected as sources in subsequent
  photometry iterations.  They are then often removed by the vetting process.
  This can bias the PSF fluxes of bright stars because the residuals identified
  as separate sources should instead be incorporated into the flux of the
  original source; see :doc:`issues`.

- **The source tables differ slightly between runs** on the same images; see
  :doc:`issues`.

- **A later step reports that a source table "is out of date"**. The image was
  reduced (or solved) again after it was measured, meaning that you should rerun
  the photometry.

Using the pipeline in Python
----------------------------

.. code-block:: python

   from pipnick.pipelines.photometry import photometry_all

   tables = photometry_all('n1_2x2_slow_full.ecsv')

See :func:`~pipnick.pipelines.photometry.photometry_all` for the full set of
parameters.
