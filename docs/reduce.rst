.. _reduce:

Basic Image Processing
======================

``pipnick_reduce`` turns the raw frames of one setup into reduced images: it
subtracts the overscan and the bias, divides the science frames by the flat
field, masks bad pixels, saturated pixels, and cosmic rays, and combines the
frames grouped in the metadata table into the combined calibration frames and
the reduced science images.  Every reduced image
carries the uncertainty of each pixel and a record of how it was made.

The algorithm
-------------

The products are built in dependency order: the bias, then the darks and
flats (which use the bias), then the science images (which use the bias and
flat).  For each combination group (``comb_id``) of the per-setup table:

1. **Process each frame**
   (:func:`~pipnick.reduction.processing.process_frame`):

   - mask the camera's bad pixels and the saturated pixels;
   - convert the counts to electrons with the gain of the readout;
   - subtract the overscan (the mean of the overscan pixels of each row) and
     trim the image to its data section;
   - set the uncertainty: the read noise and the error in the overscan level,
     plus the uncertainty of the bias and the Poisson noise of the counts
     (except for biases), propagated through the flat-fielding;
   - subtract the bias, and divide by the flat (for science frames: the sky
     flat of the filter, or the dome flat if there is no sky flat);
   - for science frames, mask the cosmic rays (``ccdproc``'s L.A.Cosmic,
     :func:`ccdproc.cosmicray_lacosmic`).

2. **Combine the frames** of the group with a sigma-clipped mean (outliers
   rejected at 5 sigma), propagating the uncertainties
   (:func:`~pipnick.reduction.processing.combine_frames`).  Flats are scaled
   to a common level first; science frames are scaled to the median exposure
   time of the group.  The frames are **not registered**: a group must share
   a pointing (a warning is given if their coordinates differ).

3. **Measure** the products: the bias level and read noise of the biases, the
   level of the flats, and the sky level and fraction of pixels affected by
   cosmic rays of the science images.  The bias level and read noise are
   checked against the camera's expected values; see :doc:`common`.

The darks are combined, but are not yet applied to the science frames.

**Reruns.**  A product is rebuilt only if its raw frames, its metadata table
values (object, filter, exposure time, coordinates, airmass), or the
calibrations it used have changed, or with ``--overwrite``.  So, after new
frames are added to the table (``pipnick_metadata --split --append``) or the
table is edited, simply rerun the reduction.

Inputs
------

- **A per-setup metadata table** (``-t``), written by ``pipnick_metadata
  --split``; see :doc:`metadata`.  Commented-out rows are ignored.
- **The raw frames** it lists, in the raw directory recorded in the table.

Running it
----------

.. code-block:: console

    pipnick_reduce -t n3_2x2_fast_full.ecsv

The options are:

- ``-t``, ``--raw_table``: The per-setup metadata table (required).
- ``-r``, ``--rdxdir``: The reduction directory.  By default, the directory
  saved in the table or, the first time, the parent of the raw directory.
  The directory is saved to the table, so the later steps find it.  A
  reduction directory can hold the products of **one setup** only.
- ``-o``, ``--overwrite``: Rebuild every product, even if it is up to date.
- ``--no_qa``: Do not write the QA figures.

Outputs
-------

In the reduction directory:

- **The products**, one FITS file per combination group:

  - ``Bias.fits``: the combined bias;
  - ``Dark_<exptime>s.fits``: the combined, bias-subtracted darks of each
    exposure time (rounded to 0.1 s), e.g., ``Dark_120s.fits``;
  - ``SkyFlat_<filter>.fits`` and ``DomeFlat_<filter>.fits``: the combined
    flats of each filter;
  - ``<object>_<filter>_C<comb_id>_<date>.fits``: the science images, e.g.,
    ``M11_V_C096_2026-08-20.fits``, named for the most common object name of
    the group, its filter (a prime written as ``p``, e.g., ``ip`` for
    ``i'``), its combination group, and the UT date of its first frame.

  Each file has three extensions:

  - ``PRIMARY``: the image, in electrons (``BUNIT = 'electron'``; scaled to
    the exposure time in the header);
  - ``MASK``: 1 for the pixels to ignore (bad, saturated, or affected by
    cosmic rays), 0 otherwise;
  - ``UNCERT``: the 1-sigma uncertainty of each pixel, in electrons.

  The header has the object, filter, exposure time, coordinates, and airmass
  from the metadata table (a header value that differs is kept in
  ``ORIGOBJ``, ``ORIGFILT``, ``ORIGEXP``, ``ORIGRA``, ``ORIGDEC``, or
  ``ORIGAIR``), the gain (``PROCGAIN``) and read noise (``PROCRN``) used, and
  pipnick's provenance:

  - ``PNVERS``: the pipnick version;
  - ``PNSETUP``, ``PNCOMBID``, ``PNNRAW``: the setup, the combination group,
    and the number of raw frames, which are also listed in ``HISTORY`` lines;
  - ``PNBIAS``, ``PNFLAT``: the bias and flat used (``none`` if none);
  - ``PNRAWHSH``, ``PNTABH``, ``PNBIASH``, ``PNFLATH``: the hashes of the raw
    frames, the table values, and the calibrations, which decide whether the
    product is up to date;
  - ``PNBIASLV``, ``PNRDNOIS`` (biases), ``PNLEVEL`` (flats), ``PNSKY``,
    ``PNCRFRAC`` (science images): the measurements of step 3.

- ``reduce/summary.ecsv``, also logged: one row per product, with its frame
  type, combination group, filter, number of frames, status (``built``, ``up
  to date``, or ``failed``, with the ``reason``), its measurements
  (``bias_level`` in ADU, ``read_noise`` and the nominal ``nominal_rn`` in
  electrons, the flat ``level``, the ``sky`` in electrons, and the
  cosmic-ray fraction ``cr_frac``), and the ``bias`` and ``flat`` used.
- ``reduce/qa/<product>.png``: a figure of each product (the image, and the
  histogram of its pixel values, with its measurements in the title); see
  :doc:`common`.
- ``log/pipnick_reduce_<time>.log``: the log; see :doc:`common`.

The metadata table is updated only to record the reduction directory.

Troubleshooting
---------------

The verdict at the end of the log lists anything that needs attention.

- **"There is no bias"** or **"There is no <filter> flat"**: the images are
  reduced without it (``PNBIAS`` or ``PNFLAT`` is ``none``), and their
  photometry will be biased or have spatial errors.  Add the calibration
  frames to the table (``pipnick_metadata --split --append``, or restore
  commented-out rows) and rerun.
- **The bias level or read noise is outside its expected range**: check the
  readout settings (e.g., a mislabeled readout speed) and the bias frames
  before trusting the reduction or its uncertainties.
- **"The frames of ... have different pointings"**: the frames of a
  combination group are combined without registration.  Give them different
  ``comb_id`` values, and combine them as a mosaic instead (see
  :doc:`mosaic`).
- **Products "not produced by this reduction"**: files of the same setup in
  the reduction directory that the table no longer produces (e.g., after a
  regrouping).  They are not deleted, and may be out of date; remove them if
  they are not needed.
- **A reduction directory with products from another setup** is refused:
  use a different directory for each setup (``-r``).
- **Flat-fielding residuals** (e.g., dust that is not fully removed) are a
  known limitation; see :doc:`issues`.

Using the pipeline in Python
----------------------------

.. code-block:: python

    from pipnick.pipelines.reduction import reduce_all

    science = reduce_all('n3_2x2_fast_full.ecsv', rdxdir='rdx_fast')

See :func:`~pipnick.pipelines.reduction.reduce_all`.
