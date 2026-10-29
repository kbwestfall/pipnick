.. _report:

Reports and inspection
======================

Two scripts collect what every step did with a night's data, without changing
or rebuilding anything: ``pipnick_report`` writes a browsable HTML report of
the whole data set, and ``pipnick_inspect`` prints the provenance of one image
and traces one of its sources.  Both read the same products as
``pipnick_status`` (see :doc:`status`): the headers, the source and
calibration tables, the summary tables, and the QA figures.

pipnick_report
--------------

.. code-block:: console

   pipnick_report -t n3_2x2_fast_full.ecsv

The report is written to the ``report`` subdirectory of the reduction
directory; open ``report/index.html`` in a browser.  Its pages are:

- ``index.html``: the state of every step for every science image (as in
  ``pipnick_status``), colored by state, with the images that need attention,
  the next command to run, and the logs, newest first;
- ``images/<image>.html``: the provenance of each science image, with its QA
  figures (see below);
- ``reduce.html``, ``astrometry.html``, ``photometry.html``,
  ``calibrate.html``, and ``mosaic.html``: the summary table and the QA
  figures of each step.  The calibration page lists the bands whose
  extinction was not measured, as in the banner at the end of the
  ``pipnick_calibrate`` log.

The links are relative, so the report works from a copy of the whole
reduction directory, but not on its own.  Products and QA figures that do not
exist are noted, and nothing is rebuilt: rerun the step (or remove
``--no_qa``) to add them, then rerun ``pipnick_report``.  The report is
written only by ``pipnick_report``, so it is out of date after any later run
of a step.

The options are:

- ``-t``, ``--raw_table``: The per-setup metadata file (required).
- ``-r``, ``--rdxdir``: The reduction directory.  If not provided, the
  directory saved in the metadata file is used.
- ``--mosaic_dir``: The directory with the mosaics, if they were written
  outside the ``mosaic`` subdirectory of the reduction directory (see
  ``pipnick_mosaic --outdir``).

The provenance of an image
--------------------------

The page of each image, and the output of ``pipnick_inspect``, give:

- **Reduction**: the state of the reduction, the object, filter, exposure
  time, airmass, night, setup, combination group, and raw frames; the bias
  and flat used; the sky level and cosmic-ray fraction; and the pipnick
  version.
- **Astrometry**: the state of the solution, its method and reference, the
  number of matched sources, the RMS, the pointing error, and the
  astrometry.net submission.
- **Photometry**: the state of the source table, the number of sources, the
  PSF model and its FWHM, the aperture radius, and the number of sources
  removed by the vetting, by reason.
- **Calibration**: the state of the calibration in the header, the
  calibration file, the zero point and color coefficient, how the extinction
  was set, and the number of standard-star measurements in the image (used,
  rejected, and unusable).
- **Mosaic**: the state of the mosaic, its label and WCS group, and the image's
  flux scale, magnitude offset, PSF FWHM before and after the matching,
  kernel, and weight in the mosaic.
- **Attention**: the reasons the image needs attention, as in
  ``pipnick_status``.

and the QA figures of every step.

pipnick_inspect
---------------

.. code-block:: console

   pipnick_inspect -t n3_2x2_fast_full.ecsv M11_V_C096_2026-08-20
   pipnick_inspect -t n3_2x2_fast_full.ecsv 96 --id 12 --display
   pipnick_inspect -t n3_2x2_fast_full.ecsv 96 --source 282.77,-6.27
   pipnick_inspect -t n3_2x2_fast_full.ecsv 96 --source 18:51:04.8,-06:16:12 --ofile c96.png

``pipnick_inspect`` prints the provenance of one image, given by its file
name (with or without the extension) or its combination group.

With ``--id`` (the ``id`` in the source table of the image) or ``--source``
(the RA and Dec, in degrees or sexagesimal with the RA in hours, matched
within ``--radius`` arcsec, default 1), it also traces one source:

- in the source table: its position, fluxes, instrumental magnitude, and
  flags;
- if no source is found at the coordinates, in the sources removed by the
  photometry (``photometry/<image>_removed.ecsv``), with the reason it was
  removed;
- in the calibrated source table: its calibrated magnitude and color, and its
  magnitudes in the other images of the same object;
- for a standard star, in the standard-star measurements of the night: the
  star, and whether the measurement was used, rejected as an outlier, or
  unusable, and why.

A source can only be found by its coordinates in an image with an astrometric
solution.

``--display`` shows the image in an interactive window, with its sources
labeled by their ``id`` (so that one can be traced with ``--id``), the removed
sources (red crosses), and the traced source (yellow circle).  ``--ofile``
writes the same figure to a file instead, which works without a display.
