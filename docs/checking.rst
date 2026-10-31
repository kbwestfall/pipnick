.. _checking:

Checking the results: pipnick_status, pipnick_report, pipnick_inspect
=====================================================================

Three scripts collect what every step did with a night's data, without
changing or rebuilding anything:

- ``pipnick_status`` prints the state of every step for every science image,
  and the next command to run;
- ``pipnick_report`` writes a browsable HTML report of the whole data set,
  with every summary and QA figure;
- ``pipnick_inspect`` prints the provenance of one image, and traces one of
  its sources through the products.

They read the products of the processing steps: the metadata table, the
headers of the reduced images, the source and calibration tables, the
summaries, and the QA figures.

pipnick_status
--------------

``pipnick_status`` answers the question "where am I, and what's broken?"

.. code-block:: console

    pipnick_status n3_2x2_fast_full.ecsv

**Inputs.**  A per-setup metadata table and its reduction directory, at any
stage of the processing.

**Options.**

- ``raw_table``: The per-setup metadata table (required).
- ``-r``, ``--rdxdir``: The reduction directory.  If not provided, the
  directory saved in the table is used.
- ``--mosaic_dir``: The directory with the mosaics, if they were written
  outside the ``mosaic`` subdirectory of the reduction directory
  (``pipnick_mosaic --outdir``).
- ``--ofile``: Also write the status table to this ECSV file.

**Output.**  A table, logged, with one row per science image: its object,
filter, exposure time, and airmass; the state of the reduction
(``reduced``), astrometric solution (``wcs``, with its RMS ``rms``),
photometry (``phot``, with the number of sources ``nsrc`` and the FWHM
``fwhm``), photometric calibration (``calibrated``, with the zero point
``zp``), and mosaic (``mosaic``); and the reasons it needs attention.  Each
step is:

- ``ok``: up to date;
- ``stale``: its inputs have changed, so it must be run again;
- ``missing``: not yet run;
- ``failed``: the step failed for this image;
- ``n/a``: the step does not apply (e.g., the image is not in a mosaic).

The ``attention`` column flags images reduced by another pipnick version
(``version``), without a bias or flat (``no bias``, ``no flat``), or with an
FWHM or astrometric RMS outside its expected range (``fwhm``, ``wcs rms``;
see :ref:`common`).  The status ends with the next command to run.

The status uses the same checks as the steps, except for the mosaics: it can
check that an image is unchanged since its mosaic was made, but not the
settings of the mosaic, which only ``pipnick_mosaic`` can check.

pipnick_report
--------------

.. code-block:: console

    pipnick_report n3_2x2_fast_full.ecsv

**Inputs.**  A per-setup metadata table and its reduction directory.

**Options.**

- ``raw_table``: The per-setup metadata table (required).
- ``-r``, ``--rdxdir``: The reduction directory.  If not provided, the
  directory saved in the table is used.
- ``--mosaic_dir``: As for ``pipnick_status``.

**Outputs.**  The report is written to the ``report`` subdirectory of the
reduction directory; open ``report/index.html`` in a browser.  Its pages are:

- ``index.html``: the status of every step for every science image (as in
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

The provenance of an image
--------------------------

The page of each image in the report, and the output of ``pipnick_inspect``,
give:

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

    pipnick_inspect n3_2x2_fast_full.ecsv M11_V_C096_2026-08-20
    pipnick_inspect n3_2x2_fast_full.ecsv 96 --id 12 --display
    pipnick_inspect n3_2x2_fast_full.ecsv 96 --source 282.77,-6.27
    pipnick_inspect n3_2x2_fast_full.ecsv 96 --source 18:51:04.8,-06:16:12 --ofile c96.png

**Inputs.**  A per-setup metadata table, its reduction directory, and one of
its science images, given by its file name (with or without the extension) or
its combination group.

**Options.**

- ``raw_table``, ``-r``, ``--rdxdir``, ``--mosaic_dir``: As for
  ``pipnick_status``.
- ``--id``: Trace the source with this ``id`` in the source table of the
  image.
- ``--source``: Trace the source at these coordinates, as RA,Dec in degrees,
  or sexagesimal with the RA in hours.  The image must have an astrometric
  solution.
- ``--radius``: The match radius in arcsec for ``--source`` (default 1).
- ``--display``: Show the image with its sources in an interactive window.
- ``--ofile``: Write the same figure to a file instead (this works without a
  display).

**Outputs.**  The provenance of the image (see above), logged.  With ``--id``
or ``--source``, the trace of the source:

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

With ``--display`` or ``--ofile``, a figure of the image with its sources
labeled by their ``id`` (so that one can be traced with ``--id``), the removed
sources (red crosses), and the traced source (yellow circle).
