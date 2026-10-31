.. _issues:

Known issues
============

The following provides a set of known issues for the current version of pipnick.
Please `submit an issue <https://github.com/UCObservatories/pipnick/issues>`__
if you find additional items!

Reduction
---------

- **Flat-fielding residuals.**  The reduction does not always remove all of
  the flat-field structure: a dust particle can remain, or the sensitivity can
  still vary slightly across the field.

Metadata
--------

- **The C2 fast readout.**  The header value (``READ-SPD``) of the fast readout
  of the Nickel C2 camera is not yet confirmed, so frames taken with it may be
  rejected by ``pipnick_metadata`` as having an undefined readout setting.
  Restore their rows in the table only if you know the gain and read noise of
  the readout are those of the slow readout; otherwise, the uncertainties of the
  reduced images are wrong.

Astrometry
----------

- **Sparse fields.**  Images with very few stars may not be solved; see the
  troubleshooting of :ref:`astrometry`.

Photometry
----------

- **Residual detections around bright stars.**  The residuals of the PSF
  model around bright stars can be detected as sources.  These are removed by
  the vetting, but the fits of the PSF magnitudes of nearby sources can be
  biased by them.  Compare the PSF and aperture magnitudes of sources near
  bright stars.

- **Run-to-run differences.**  The source tables of the same image can differ
  slightly between runs of ``pipnick_photometry`` (e.g., in the faintest
  detections).  Rebuilding the photometry of an image, even with the same
  settings, can therefore change the calibration and mosaics that use it.

- **Removed sources that appear to be outside the image.**  Some sources
  removed as ``failed`` have starting coordinates (``x_init``, ``y_init`` in
  ``photometry/<image>_removed.ecsv``) far outside the image.  They were
  detected within the image, but their fit in that iteration ran off the image.  With ``--mode all``, the next
  iteration starts their fit from that position, where there are no pixels to
  fit, so the fit fails.  These sources do not affect the source table.

Mosaics
-------

- **One unsolved image stops every mosaic.**  ``pipnick_mosaic`` requires an
  astrometric solution and a PSF model for every image assigned to a mosaic, and
  stops with an error otherwise, instead of skipping the affected mosaic.  For
  now, you must clear the labels of the problematic image (see :ref:`mosaic`).

Status and report
-----------------

- **Bands that cannot be calibrated.**  ``pipnick_status`` and
  ``pipnick_report`` suggest running ``pipnick_calibrate`` for an image whose
  band has no calibration, even when the band cannot be calibrated (e.g.,
  there are no standard stars in the band, or the filter has no standard
  magnitudes).  The ``pipnick_calibrate`` log says why the band was not
  calibrated; if nothing has changed, rerunning it does not help.
