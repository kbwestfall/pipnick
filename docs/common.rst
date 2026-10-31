.. _common:

Common script behavior
======================

The pipnick processing scripts share a number of commonalities.  Discussions of
these are collected here to minimize repetition in the documentation.

Command-line options
--------------------

All pipnick scripts will provide basic help when you use the ``-h`` option, such
as:

.. code-block:: bash

    pipnick_metadata -h

This lists all the required and optional arguments and a short description of
their purpose.


The reduction directory
-----------------------

Reductions and analysis are expected to be performed for each mode of the
instrument for each night in its own top-level reduction directory.  For
simplicity, pipnick adopts a fixed directory structure below this top-level
directory.

.. warning::

    If you move any files within the following directory structure, pipnick will
    assume it was not produced.  Execution of subsequent processing scripts may
    try to reproduce moved files (adding them to their expected on-disk
    locations), or they may not be able to proceed because of the missing file.

When you run ``pipnick_reduce``, the reduction directory is saved in the header
of the ECSV metadata table, and it is referenced this way through the other
processing scripts.

The directory structure pipnick uses for its output is as follows:

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

Completion of each processing step can result in both new files and edits to the
reduced science images.  For example, ``pipnick_astrometry`` not only adds the
new data in the ``astrometry/`` directory, but also adds the astrometric solution
(WCS) to the processed file headers.

Logs
----

Each processing script writes a log file by default.  For example, the log file
for the basic image processing is saved to
``<rdxdir>/log/pipnick_reduce_<UT-timestamp>.log`` (e.g.,
``pipnick_reduce_20261028T153614.log``).  A new log file is produced every
time the script is executed, meaning logs are never overwritten or appended to.
If two logs would have the same name, the second one will be named
``..._2.log``, etc.

Scripts that do not operate within the reduction directory will write their logs
in a ``log/`` subdirectory nearest to their output: ``pipnick_metadata`` writes
logs to ``log/`` in the default reduction directory (the parent of the raw
directory), ``pipnick_mosaic`` writes logs to ``log/`` in the ``--outdir``, if
it is set, and ``pipnick_rgb`` writes logs to ``log/`` in the directory of its
output image.

Every processing script has the same logging options:

- ``-v``, ``--verbosity``: This sets the verbosity of the output written to the
  terminal: 0 only prints warnings and errors, 1 adds informational messages
  (the default), and 2 adds debugging messages and the python calling sequence.

- ``--log_file``: Use this to define the log file name directly instead of using
  the default naming convention.  **The log file must not exist.**

- ``--log_level``: Use this to set the verbosity of the log file if you want it
  to be different from the terminal verbosity (the default).

Versions
--------

The first line of every log is the full pipnick version.  It is possible for
different pipnick versions to be used for different steps.  For this reason, the
version is saved in multiple places with names specific to the step executed:

- the ``PNVERS`` keyword is added to the header of
  the reduced images and mosaics;

- the ``PNASTVER`` keyword is the pipnick version specific to the astrometric
  solution, and the ``PNCALVER`` keyword is specific to the photometric
  calibration;

- the ``pipnick_version`` is added to the metadata tables, source tables,
  calibration files, and summaries.

When a step reuses an existing product made by a different version of pipnick,
it issues one warning listing these products.  They are still used: the version
never causes a product to be rebuilt.  **However, we strongly recommend
rerunning all of the processing steps anytime you upgrade pipnick.**  This can
be done using the ``--overwrite`` option (``--resolve`` for
``pipnick_astrometry``), by deleting existing data, or by
sending the output to a new reduction directory.

.. note::
    
    The ``pipnick_version`` script always reports the currently installed version.

Outcome status and summaries
----------------------------

Each step ends with an outcome status indicating the end state of the execution:

- ``<STEP> OK``: Everything executed successfully

- ``<STEP> NEEDS ATTENTION``: The processing/analysis step finished, but some
  results need to be inspected.  This could be a missing flat or a measurement
  outside its expected range.  The specific problems encountered are repeated in
  the outcome.

- ``<STEP> FAILED``: A fatal error occurred, meaning that outputs were not fully
  constructed.  The specific failures are listed and the script ends with a
  non-zero exit status.

The outcome status is followed by the next command to run.  Each step also
writes a summary table to its output directory, ``<step>/summary.ecsv``, with
one row per product and its status; the summary is also logged.  The columns of
each summary are described on the page of its script.

Reruns
------

Every step can be rerun at any time, and only does the work that is needed:

- a product is rebuilt only if its inputs or settings changed (e.g., new raw
  frames, an edited metadata table, a new astrometric solution, other
  options), which is recorded by hashes in its header or metadata;
- a product that is up to date is reported as such, and reused;
- ``--overwrite`` (``--resolve`` for ``pipnick_astrometry``) rebuilds every
  product of the step.

After adding frames to a night or changing a setting, you can simply rerun the
steps in order, and each step will only rebuild the outputs that must be
updated.  The ``pipnick_status`` script will report which steps are out of date
(see :ref:`checking`).

QA figures
----------

Each step writes a set of fixed-format QA figures to the ``qa`` subdirectory of
its main directory, such as ``reduce/``.  You can skip writing these figures
using the ``--no_qa`` option; however, you are strongly encouraged to keep them so that
you can check the results.

- ``reduce/qa/<product>.png``: An image of each processed frame and a histogram
  of its pixel values.

- ``astrometry/qa/<image>.png``: The processed image with the detected sources
  overlaid.  For each image solved, residuals between the source and solution
  coordinates are shown.

- ``photometry/qa/<image>_psf.png`` and ``photometry/qa/<image>_sources.png``:
  Quality check for the PSF model and the measured sources.

- ``calibrate/qa/<setup>_<night>_residuals.png``: A plot of the calibration
  residuals.

- ``mosaic/qa/<mosaic>.png``: An image of the mosaic, its sky coverage, the flux
  scales, and the PSF matching.

Each QA figure is redrawn if its output product is rebuilt, and it is drawn for
any product that lacks one.  A figure that cannot be drawn is flagged as needing
attention, but does not stop the step.  ``pipnick_report`` collects all of the
figures in one place (see :ref:`checking`).

Expected ranges
---------------

Specific measurements are checked against the expected values for the camera, as
set in its configuration file (``[expected]`` and ``[quality]``):

- the bias level (ADU);

- the read noise, against the nominal value of the readout;

- the PSF FWHM of each image (arcsec);

- the RMS of each astrometric solution (arcsec).

A value outside its range is flagged as needing attention, but data are never
rejected due to these checks.
