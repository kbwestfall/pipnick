.. _common:

Common behaviour of the scripts
===============================

The processing scripts share their layout of the outputs, their logs and
verdicts, their rules for rerunning, and their QA figures.  These are
described here once, and referred to by the page of each script.

The reduction directory
-----------------------

Each per-setup metadata table is reduced into its own **reduction
directory**, which can hold the products of only one setup.  It is set by
``pipnick_reduce -r`` (by default, the parent of the raw directory), saved to
the table, and found there by the later steps.  Its layout is:

.. code-block:: text

    <rdxdir>/
        Bias.fits, DomeFlat_V.fits, ...        the calibration frames (pipnick_reduce)
        M11_V_C096_2026-08-20.fits, ...        the reduced science images
        reduce/                                summary.ecsv, qa/
        astrometry/                            the astrometry.net and Gaia caches,
                                               summary.ecsv, qa/, dry_run/
        photometry/                            <image>_phot.ecsv, <image>_removed.ecsv,
                                               summary.ecsv, qa/
        calibrate/                             <setup>_<night>.toml, <image>_cal.ecsv,
                                               <setup>_<night>_residuals.ecsv,
                                               summary.ecsv, qa/
        mosaic/                                <mosaic>.fits, <wcs_group>_wcs.fits,
                                               summary.ecsv, qa/
        report/                                index.html, ... (pipnick_report)
        log/                                   the logs of every script

The reduced science images are the hub: each later step adds to their headers
(the astrometric solution, the photometric calibration) or writes products
named for them.

Logs
----

Each script writes a log by default, to a new file in the ``log``
subdirectory of the reduction directory, named for the script and the time
(UT) it started; e.g., ``<rdxdir>/log/pipnick_reduce_20261028T175048.log``.
Logs are never overwritten or appended to: two runs started in the same
second get distinct names (``..._2.log``).

Scripts without a reduction directory write their logs to the nearest output
directory: ``pipnick_metadata`` to the default reduction directory (the
parent of the raw directory), ``pipnick_mosaic`` with ``--outdir`` to that
directory, and ``pipnick_rgb`` next to its output image.

Every script has the same logging options:

- ``-v``, ``--verbosity``: The verbosity of the output to the terminal: 0
  (warnings and errors), 1 (the default; also informational messages), or 2
  (also debugging messages and the calling sequence).
- ``--log_file``: A different log file, which must not exist.
- ``--log_level``: The verbosity of the log file, if it differs from the
  terminal's.

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
``--overwrite`` is recommended.  ``pipnick_version`` reports the installed
version.

Verdicts and summaries
----------------------

Each step ends with a verdict:

- ``<STEP> OK``: everything worked;
- ``<STEP> NEEDS ATTENTION``: the step finished, but some results need a
  look, e.g., a missing flat or a measurement outside its expected range;
  the problems are repeated in the verdict;
- ``<STEP> FAILED``: some products could not be constructed.  The other
  products are still constructed, the failures are listed, and the script
  ends with a non-zero exit status.

The verdict is followed by the next command to run.  Each step also writes a
summary table to its output directory, ``<step>/summary.ecsv``, with one row
per product and its status; the summary is also logged.  The columns of each
summary are described on the page of its script.

Reruns
------

Every step can be rerun at any time, and only does the work that is needed:

- a product is rebuilt only if its inputs or settings changed (e.g., new raw
  frames, an edited metadata table, a new astrometric solution, other
  options), which is recorded by hashes in its header or metadata;
- a product that is up to date is reported as such, and reused;
- ``--overwrite`` rebuilds every product of the step.

So, after adding frames to a night or changing a setting, rerun the steps in
order; each one rebuilds only what changed.  ``pipnick_status`` reports which
steps are out of date (see :doc:`checking`).

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
that lacks one.  A figure that cannot be drawn needs attention, but does not
stop the step.  ``pipnick_report`` collects all of the figures in one place
(see :doc:`checking`).

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
