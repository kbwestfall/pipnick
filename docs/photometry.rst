Point-Source and Aperture Photometry
====================================

``pipnick_photometry`` measures the sources in each reduced science image with
point-spread-function (PSF) and aperture photometry, and writes a source table
for each image.  It does not require an astrometric solution: the pixel
coordinates of the sources are always provided, and their celestial
coordinates are provided for images with a solution from
``pipnick_astrometry`` (see :doc:`astrometry`).

Overview
--------

For each image:

1. **Model the PSF.**  Bright, isolated stars are detected, background-
   subtracted cutouts around them are stacked (``--stack``, mean by default),
   and a circular or elliptical (``--fittype``) Moffat function is fit to the
   stack.  The fit sets the FWHM used for all the scales below.

2. **Subtract the background.**  A smooth, 2D background is subtracted
   from the image.

3. **PSF photometry.**  Sources are detected above ``--thresh`` times the
   background noise (with ``photutils``' ``DAOStarFinder``), no closer than
   ``--min_sep`` FWHM to each other, and their fluxes and positions are fit
   with the PSF, its shape fixed.  Sources closer than ``--group_sep`` FWHM are
   fit together, and detection and fitting are repeated on the residual image
   (``photutils``' ``IterativePSFPhotometry``; see ``--mode``).  With
   ``--group``, sources fit together are consolidated into one.

4. **Vet the sources.**  Sources whose fits failed, whose fitted positions
   are outside the image or moved more than one FWHM from where they were
   detected, or whose fitted fluxes are not positive are removed.  These are
   typically spurious detections in the residuals of bright stars.  They are
   listed, with the reason, in ``photometry/<image>_removed.ecsv``.

5. **Aperture photometry.**  The background-subtracted flux in a circular
   aperture of radius ``--aper_radius`` FWHM is measured at the position of
   each source.

The pixel errors used to weight the fits and compute the flux errors are the
uncertainties propagated by ``pipnick_reduce`` (the ``UNCERT`` extension of
the reduced image): the read noise (including the error in the subtracted
overscan level), the Poisson noise of the counts, and the uncertainties of the
bias and the flat, propagated through the combination of frames.  The object,
filter, exposure time, and airmass are taken from the metadata table.

Using the command-line script
-----------------------------

.. code-block:: console

   pipnick_photometry -t n1_2x2_slow_full.ecsv

The options are:

- ``-t``, ``--raw_table``: The per-setup metadata file used to reduce the
  images (required).
- ``-r``, ``--rdxdir``: The reduction directory.  If not provided, the
  directory saved in the metadata file is used.
- ``--images``: Only measure these images, given by file name (with or without
  the extension) or combination group (``comb_id``).
- ``-o``, ``--overwrite``: Measure the images again, even if they are up to
  date.
- ``--no_qa``: Do not write the QA figures to ``photometry/qa``: the PSF model
  (``<image>_psf.png``: the stacked stars, the model, the residuals, and the
  radial profile) and the measured sources (``<image>_sources.png``).  By
  default, the figures are redrawn for each measured image, and drawn for an
  image that is up to date but lacks them (its PSF stars are stacked again
  from the image, without measuring the sources again).
- ``--thresh``: The source detection threshold, as a multiple of the
  background noise (default 8).
- ``--psf_thresh``: The detection threshold for the stars used to model the
  PSF (default 15).
- ``--fwhm``: The initial guess for the FWHM in arcsec (default: typical for
  the camera).
- ``--fittype``: The PSF model, ``circ`` (default) or ``ellip``.
- ``--stack``: How the PSF stars are stacked, ``mean`` (default) or
  ``median``.
- ``--mode``: The iterative fitting mode, ``all`` (default) or ``new``.
- ``--min_sep``, ``--group_sep``: The minimum separation of detected sources,
  and the separation below which sources are fit together, in FWHM (defaults
  1 and 2).
- ``--group``: Consolidate the sources fit together into single sources.
- ``--aper_radius``: The radius of the photometric apertures in FWHM (default
  2).

Rerunning the script is cheap: an image is measured again only if it changed
(e.g., it was reduced again) or the settings changed.  If only its
astrometric solution changed, for example because ``pipnick_astrometry`` was
run after ``pipnick_photometry``, only the celestial coordinates are updated.

Outputs
-------

The source table for each image is ``photometry/<image>_phot.ecsv`` in the
reduction directory, with these columns:

- ``id``, ``group_id``, ``group_size``, ``iter_detected``: The source, the
  group of sources fit together, and the iteration in which it was detected.
- ``x``, ``y``, ``x_err``, ``y_err``: The 0-based pixel coordinates.  These
  are always provided.
- ``ra``, ``dec``: The ICRS coordinates in degrees, if the image has an
  astrometric solution from ``pipnick_astrometry``, and masked otherwise.  The
  nominal WCS that the reduced images inherit from the raw headers is never
  used.
- ``flux_psf``, ``flux_psf_err``: The PSF flux in electrons.
- ``flux_aper``, ``flux_aper_err``: The aperture flux in electrons, masked if
  the aperture includes masked pixels or extends beyond the image.
- ``mag_inst``, ``mag_inst_err``: The instrumental magnitude,
  :math:`-2.5 \log_{10}(f / t)` for the PSF flux :math:`f` and the exposure
  time :math:`t`.
- ``local_bkg``, ``qfit``, ``cfit``, ``reduced_chi2``, ``flags``: The local
  background and the fit-quality measures and flags of ``photutils``.

Values that cannot be measured are masked.  The table metadata record the
image and its object, filter, exposure time, airmass, and MJD; the PSF model
(type, parameters, FWHM in pixels and arcsec, and the number of stars used);
the number of sources removed by the vetting; the settings; and the
provenance of the image and its astrometric solution.

The sources removed by the vetting are written to
``photometry/<image>_removed.ecsv``, with the 0-based pixel coordinates where
they were detected (``x_init``, ``y_init``), their fitted coordinates and flux
(``x_fit``, ``y_fit``, ``flux_fit``; masked if the fit failed), the iteration
of the PSF photometry that detected them (``iter_detected``), and the reason
they were removed (``failed``, ``outside``, ``moved``, or ``nonpositive``).
Source tables measured before this file existed gain it when they are next
measured.

``photometry/summary.ecsv``, also logged, summarizes each image: its status
(``measured``, ``up to date``, ``coordinates updated``, ``missing``, or
``failed``), the number of sources, the FWHM, the number of PSF stars, and the
astrometric solution method (``none`` if the coordinates are masked).

The instrumental magnitudes are not calibrated; ``pipnick_calibrate`` fits the
photometric calibration from the Landolt standard stars (see :doc:`calibrate`).

Using the pipeline in Python
----------------------------

.. code-block:: python

   from pipnick.pipelines.photometry import photometry_all

   tables = photometry_all('n1_2x2_slow_full.ecsv')

See :func:`~pipnick.pipelines.photometry.photometry_all` for the full set of
parameters.
