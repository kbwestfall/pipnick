.. _status:

Logs, verdicts, and status
==========================

Every pipnick script keeps a record of what it did, ends with a verdict on
the result, and says what to run next.  ``pipnick_status`` gives the state of
every step for every science image of a night.

Logs
----

Each script writes a log by default, to a new file in the ``log``
subdirectory of the reduction directory, named for the script and the time
(UT) it started; e.g., ``<rdxdir>/log/pipnick_reduce_20261028T175048.log``.
Logs are never overwritten or appended to: two runs started in the same
second get distinct names (``..._2.log``).  Use ``--log_file`` to choose a
different file, which must not exist.

Scripts without a reduction directory write their logs to the nearest output
directory: ``pipnick_metadata`` to the default reduction directory (the
parent of the raw directory), ``pipnick_mosaic`` with ``--outdir`` to that
directory, and ``pipnick_rgb`` next to its output image.

Versions
--------

The first line of every log is the full pipnick version, and the version is
recorded in every product:

- the ``PNVERS`` header keyword of the reduced images and mosaics;
- the ``PNASTVER`` keyword of an astrometric solution, and the ``PNCALVER``
  keyword of a photometric calibration recorded in a header;
- the ``pipnick_version`` metadata of the metadata tables, source tables,
  calibration files, and summaries.

When a step reuses an existing product made by a different version of pipnick
(or by an earlier version that did not record it), it issues one warning
listing these products.  They are still used: the version never causes a
product to be rebuilt.  However, rerunning all of the processing steps with
``--overwrite`` is recommended.

Verdicts
--------

Each step ends with a verdict:

- ``<STEP> OK``: everything worked;
- ``<STEP> NEEDS ATTENTION``: the step finished, but some results need a
  look, e.g., a missing flat or a measurement outside its expected range;
  the problems are repeated in the verdict;
- ``<STEP> FAILED``: some products could not be constructed.  The other
  products are still constructed, the failures are listed, and the script
  ends with a non-zero exit status.

The verdict is followed by the next command to run.  Each step also writes a
summary table to its output directory, e.g., ``reduce/summary.ecsv`` for the
reduction, which gives the bias level and read noise of the bias, the levels
of the flats, and the sky level and cosmic-ray fraction of each science image.

QA figures
----------

Each step writes QA figures by default, to the ``qa`` subdirectory of its
output directory (use ``--no_qa`` to skip them):

- ``reduce/qa/<product>.png``: each combined calibration frame and reduced
  science image, with the histogram of its pixel values and its measurements;
- ``astrometry/qa/<image>.png``: the detected sources and, for an image solved
  in the run, the sources matched by the solution with their residuals;
- ``photometry/qa/<image>_psf.png`` and ``<image>_sources.png``: the PSF model
  and the measured sources;
- ``calibrate/qa/<setup>_<night>_residuals.png``: the residuals of the
  calibration;
- ``mosaic/qa/<mosaic>.png``: the mosaic, its coverage, the flux scales, and
  the PSF matching.

A figure is redrawn when its product is rebuilt, and drawn for any product
that lacks one, e.g., a product that is up to date but was made before the
figures were written by default.  A figure that cannot be drawn needs
attention, but does not stop the step.

Expected ranges
---------------

These measurements are checked against the expected values for the camera,
set in its configuration file (``[expected]`` and ``[quality]``):

- the bias level (ADU);
- the read noise, against the nominal value of the readout;
- the PSF FWHM of each image (arcsec);
- the RMS of each astrometric solution (arcsec).

A value outside its range needs attention; data are never rejected because of
it.

pipnick_status
--------------

``pipnick_status -t <table>`` reports, for each science image of a per-setup
metadata table, the state of the reduction, astrometry, photometry,
calibration, and mosaic, with the astrometric RMS, the number of sources, the
FWHM, and the zero point.  Each step is:

- ``ok``: up to date;
- ``stale``: its inputs have changed, so it must be run again;
- ``missing``: not yet run;
- ``failed``: the step failed for this image;
- ``n/a``: the step does not apply (e.g., the image is not in a mosaic).

The ``attention`` column flags images reduced by another pipnick version
(``version``), without a bias or flat (``no bias``, ``no flat``), or with an
FWHM or astrometric RMS outside its expected range (``fwhm``, ``wcs rms``).
The status ends with the next command to run.  It changes nothing; use
``--ofile`` to also write the table to an ECSV file.
