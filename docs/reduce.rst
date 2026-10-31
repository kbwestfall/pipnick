.. _reduce:

Basic Image Processing
======================

Using the raw data listed and organized in the metadata table,
``pipnick_reduce`` performs basic image processing of the raw data, identifying
and removing instrumental effects.  Every reduced image carries the uncertainty
of each pixel and a record of how it was made.

The algorithm
-------------

Image processing produces the calibrations first and then processes the science
frames.  The algorithm is as follows:

1. **Process each individual frame**
   (:func:`~pipnick.reduction.processing.process_frame`):

   - Mask the camera's bad pixels and the saturated pixels.

   - Convert the counts to electrons using the hard-coded detector gain.

   - Subtract the overscan, as determined using the mean of the overscan pixels
     of each row, and trim the image to its data section.

   - Subtract the bias.
   
   - Divide by the flat.  Sky flats are used if available; dome flats are used
     otherwise.

   - For science frames, identify and mask cosmic rays, using
     :func:`ccdproc.cosmicray_lacosmic`.

   - Estimate the uncertainty starting with the hard-coded detector readnoise
     and propagating errors through the other processing steps.

2. **Combine frames** (:func:`~pipnick.reduction.processing.combine_frames`):
   Frames that have the same ``comb_id`` value are combined using a
   sigma-clipped mean (5-sigma outliers are rejected); uncertainties are
   propagated.  Flats are scaled to a common level before combining.  Science
   frames should nominally all have the same exposure time; when this is not
   true, the images are scaled to the median exposure time of the group before
   combining.  The frames are **not registered spatially**: this step combines
   images pixel-by-pixel.

3. **Measure quantitative metrics**: As feedback for the image processing, the
   code reports the bias level and read noise of the biases, the level of the
   flats, and the sky level and fraction of pixels affected by cosmic rays of
   the science images.  The bias level and read noise are checked against the
   camera's expected values; see :ref:`common`.

Dark images are combined, but they are not yet included in the image processing.

If ``pipnick_reduce`` is executed more than once, processing is only redone its
raw frames, its metadata table values (object, filter, exposure time,
coordinates, airmass), or the calibrations it used have changed.  It is possible
to force all the data to be reprocessed using ``--overwrite``.

Inputs
------

The primary inputs for ``pipnick_reduce`` are:

- **The ECSV metadata table** (``raw_table``), written by ``pipnick_metadata``;
  see :ref:`metadata`.

- **The raw frames**

Execution
---------

An example call looks like this:

.. code-block:: console

    pipnick_reduce n3_2x2_fast_full.ecsv

The required argument is the per-setup metadata table.

The options are:

- ``-r``, ``--rdxdir``: The path for the reduction output files.  If not
  provided and not set in the metadata table, the parent directory of the
  directory containing the raw data is used.  We **strongly** recommend you set
  this directory explicitly.  The directory is saved to metadata table so that
  it is not a required command-line argument for all remaining scripts.

- ``-o``, ``--overwrite``: Rebuild every product, even if it is up to date.

- ``--no_qa``: Do not write the QA figures.

.. _reduce_outputs:

Outputs
-------

The reduction directory will contain the following processed FITS files:

- ``Bias.fits``: The combined bias.

- ``Dark_<exptime>s.fits``: The combined, bias-subtracted darks of each
    exposure time (rounded to 0.1 s), e.g., ``Dark_120s.fits``.

- ``SkyFlat_<filter>.fits`` and ``DomeFlat_<filter>.fits``: The combined
    flats for each filter.

- ``<object>_<filter>_C<comb_id>_<date>.fits``:  The science images, e.g.,
    ``M11_V_C096_2026-08-20.fits``, named for the most common object name of
    the group, its filter (a prime written as ``p``, e.g., ``ip`` for
    ``i'``), its combination group, and the UT date of its first frame.

Each FITS file has three extensions:

- ``PRIMARY``:  The processed image data in electrons.

- ``MASK``: 1 for the pixels to ignore (bad, saturated, or affected by cosmic
  rays), 0 otherwise;

- ``UNCERT``: the 1-sigma uncertainty of each pixel in electrons.

The primary header contains a copy of the header from the first frame in the
combination group.  Additionally, it contains the object, filter, exposure time,
coordinates, and airmass from the metadata table; if any of these values differ
from the one in the raw header, an ``ORIG*`` version is kept as a record of the
difference.  For example, ``ORIGOBJ`` may be added to include the original
``OBJECT``; other such keywords are ``ORIGFILT``, ``ORIGEXP``, ``ORIGRA``,
``ORIGDEC``, and ``ORIGAIR``.

Additional header keywords included are:

- ``PROCGAIN``, ``PROCRN``: gain and readnoise

- ``PNVERS``: the pipnick version

- ``PNSETUP``, ``PNCOMBID``, ``PNNRAW``: the setup, the combination group,
  and the number of raw frames, which are also listed in ``HISTORY`` lines

- ``PNBIAS``, ``PNFLAT``: the bias and flat used (``none`` if these corrections
  were not applied)

- ``PNRAWHSH``, ``PNTABH``, ``PNBIASH``, ``PNFLATH``: the hashes of the raw
  frames, the table values, and the calibrations, which are used to decide
  whether the product is up to date 

- ``PNBIASLV``, ``PNRDNOIS`` (biases), ``PNLEVEL`` (flats), ``PNSKY``,
  ``PNCRFRAC`` (science images): the measurements from step 3 of the algorithm.

Additionally, the script writes a summary of the results to
``reduce/summary.ecsv``.  This is a table that includes one row per product,
with its frame type, combination group, filter, number of frames, status
(``built``, ``up to date``, or ``failed``, with the ``reason``), its
measurements (``bias_level`` in ADU, ``read_noise`` and the nominal
``nominal_rn`` in electrons, the flat ``level``, the ``sky`` in electrons, and
the cosmic-ray fraction ``cr_frac``), and the ``bias`` and ``flat`` used.

QA plots are added to ``reduce/qa`` that show a plot of the image and a
histogram of its pixel values; see :ref:`common`.

Finally, the log is saved to ``log/pipnick_reduce_<time>.log``; see
:ref:`common`.

The metadata table is updated only to record the reduction directory.

Troubleshooting
---------------

The outcome of the script is provided at the end of the log and lists anything
that needs attention.

- **"There is no bias"** or **"There is no <filter> flat"**:  There are missing
  calibration images.  The reduction will continue without the relevant
  corrections, meaning further processing will be severely limited.  We recomend
  addint the relevant raw calibration frames to the metadata table (e.g.,
  ``pipnick_metadata --split --append`` or restore commented-out rows) and
  rerun.

- **The bias level or read noise is outside its expected range**:  Check the
  readout settings (e.g., a mislabeled readout speed) and the bias frames before
  trusting the reduction or its uncertainties.

- **"The frames of ... have different pointings"**:  The frames of a combination
  group are combined without registration.  The code will use the celestial
  coordinates provided by the metadata to determine if there is a significant
  difference between the pointings of any of the frames being combined.  Check
  that the frames have identical pointings, then either correct the RA/Dec
  coordinates in the metadata table or give them different ``comb_id`` values
  and combine them as a mosaic instead (see :ref:`mosaic`).

- **Products "not produced by this reduction"**: Files were found in the
  reduction directory that were not part of this reduction process.  This may be
  because frames were regrouped and/or out of date.  In any case, they likely
  should be removed to ensure the processing output is understood.

- **A reduction directory with products from another setup**:  The code will
  refuse to write products to a directory that already has output from a
  different setup.  Use a different directory for each setup (``-r``).

- **Flat-fielding residuals** (e.g., dust that is not fully removed) are a
  known limitation; see :ref:`issues`.

Using the pipeline in Python
----------------------------

.. code-block:: python

    from pipnick.pipelines.reduction import reduce_all

    science = reduce_all('n3_2x2_fast_full.ecsv', rdxdir='rdx_fast')

See :func:`~pipnick.pipelines.reduction.reduce_all`.
