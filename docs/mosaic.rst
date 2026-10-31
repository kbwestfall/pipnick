.. _mosaic:

Image Mosaicking
================

For observations taken with small pointing offsets (on-sky dithering), pipnick
can register them to the same astrometric coordinates (WCS) and combine them
using ``pipnick_mosaic``.  Additionally, mosaics of a single patch of sky
observed in different filters can be forced to share an identical output WCS and
shape, enabling quick construction of three-color images; see :ref:`rgb`.

Building image mosaics requires completion of both the astrometry and
photometry:  The images are registered using their astrometric solutions, and
the PSF models and source tables of the photometry are used to correctly match
images in the mosaic to the same spatial resolution and flux scale.  A typical
workflow will be

.. code-block:: console

   pipnick_metadata raw n3 --split --assign_mosaics
   # Edit the n3_2x2_fast_full.ecsv file!
   pipnick_reduce n3_2x2_fast_full.ecsv
   pipnick_astrometry n3_2x2_fast_full.ecsv
   pipnick_photometry n3_2x2_fast_full.ecsv
   pipnick_mosaic n3_2x2_fast_full.ecsv

Inputs
------

- **The ECSV metadata table** (``raw_table``), written by ``pipnick_metadata``;
  see :ref:`metadata`.  The table must include the ``mosaic`` and ``wcs_group``
  columns (see below).  They can be added using ``--assign_mosaics`` or edited
  into the metadata table directly.

- **The reduced science images** (see :ref:`reduce`), with their
  **astrometric solutions** from ``pipnick_astrometry`` (see
  :ref:`astrometry`); an image without one is an error.

- **Their source tables** (see :ref:`photometry`), up to date, for the PSF
  models (unless ``--psf_match none``) and the scaling by common sources (the
  default ``--scale sources``).

- When setting ``--scale calib`` (optional), **their photometric calibrations**
  must be in the image headers (``pipnick_calibrate --apply image``; see
  :ref:`calibrate`).

Defining the mosaics
--------------------

The mosaics are defined by two columns of the per-setup metadata tables,
with *string* labels (unlike the *integer* labels used by ``comb_id``):

- ``mosaic``: the frames with the same label are combined into one mosaic;

- ``wcs_group``: the mosaics with the same label share an output WCS and shape.
  If it is empty, a mosaic is its own WCS group (i.e., no mosaics are forced to
  have the same WCS).

An empty label means the frame is not mosaicked.  The labels apply to the
reduced images: all of the frames of a combination group (``comb_id``) must
have the same labels.  A mosaic must belong to one WCS group, and a warning is
issued if it combines images with different filters.

By default, the columns are empty, i.e., no mosaics are defined.  Either edit
the tables, or assign the mosaics automatically with ``pipnick_metadata
--split --assign_mosaics``:

- The science frames are grouped by pointing: pointings within
  ``--mosaic_radius`` (by default, half of the field of view) of each other,
  directly or through other pointings, form one WCS group.

- Within each pointing group, the frames with the same filter form one mosaic.

- A filter with a single image in a pointing group is not mosaicked.

- The WCS group is labeled with the most common object name of its frames,
  without a trailing filter name, e.g. ``NGC_7331`` for ``NGC 7331 V``, and the
  mosaic with the WCS group and the filter, e.g. ``NGC_7331_V``.  The object
  names are not used to group the frames because they can be incorrect (e.g.,
  observers can mistakenly neglect to update the filter in the name).

- Only frames without labels are assigned, meaning direct edits are always kept.
  When running ``pipnick_metadata`` with ``--append``, the labels are assigned
  even if there are no new frames.

A mosaic can have a single image.  To sample several images on the same pixel
grid, e.g., one image of an object in each of three filters for a color
image, give each image its own ``mosaic`` label and all of them the same
``wcs_group``; each is then only reprojected onto the WCS of the group (it is
not convolved, unless ``--psf_match group`` is used).

The mosaic labels are global: the same label in different tables refers to the same
mosaic.  For example, a field observed with two readout modes has two setups,
and so two tables and two reduction directories; give both tables to
``pipnick_mosaic`` (with ``--outdir``) to combine its images.  When the tables
of several setups are built together, ``--assign_mosaics`` gives the same
pointing and filter the same labels in all of them.

.. note::

    For now, a mosaic can only combine images from one night (UT date).  The
    labels, the support for several tables, and the provenance recorded for each
    image are the basis for mosaics of several nights.  This functionality will
    need to be enabled explicitly such that common labels, like ``M11_V``, do
    not combine different nights by accident.

The algorithm
-------------

The output WCS of each WCS group is a north-up gnomonic projection in the ICRS
that covers all of the images of the group (in all filters), with the finest
pixel scale of the images (or ``--pixel_scale``); see
:func:`reproject.mosaicking.find_optimal_celestial_wcs`.  It is saved to
``<wcs_group>_wcs.fits`` and reused while the group is unchanged, so every
mosaic of the group has the identical WCS and shape.

Each mosaic is then constructed as follows.

1. **Prepare the images.**  Each image is converted to electrons per second, and
   its background is subtracted (``--background``).  By default the background
   is a constant as determined by the median of the background-only pixels.  A
   smooth 2D background (256-pixel boxes) can also be used; beware applying this
   to images with objects that fill much of the field-of-view.

2. **Match the PSFs** (``--psf_match``).  The seeing and focus can change during
   the observations, so by default every image is convolved to the worst seeing
   of any image in the mosaic.  The target is a circular Moffat function with
   ``--psf_margin`` (1.05) times the largest FWHM of the images, and the median
   of their Moffat power indices.  For an elliptical PSF model
   (``pipnick_photometry --fittype ellip``), the FWHM along its major axis is
   used, so that the PSF is broadened to the target in every direction.  Each
   image is convolved with the kernel that matches its PSF model (from the
   photometry) to the target, computed with
   :func:`photutils.psf_matching.make_kernel` and tapered to a circle; if that
   fails, a Gaussian kernel is used.  For constructing color images, use
   ``--psf_match group`` to match all of the mosaics of a WCS group to the worst
   seeing of the group.  ``--psf_match none`` turns the matching off.  A warning
   is issued if one image has much worse seeing than the others, because it
   degrades the resolution of the whole mosaic.

3. **Scale the fluxes** (``--scale``).  The transparency, airmass, and exposure
   time can differ between the images, so by default each image is scaled to the
   flux scale of a reference image.  By default, the reference image is the one
   observed at the lowest airmass, but it can be set directly using
   ``--reference``.  The scales are measured from the sources common to
   overlapping images: their fluxes are measured on the matched images, in
   apertures of 1.5 times the target FWHM with a local background (with
   ``--psf_match none``, their PSF fluxes from the source tables are used), and
   the magnitude offsets of all of the images are fit together.  Images that do not
   share any sources with the reference are scaled through the other images in
   the mosaic that they overlap.  Images that share too few sources with the
   others (``--min_common``) are not scaled and a warning is issued.  The
   ``--scale calib`` option will use the photometric calibrations of the images
   instead (see :ref:`calibrate`), and ``--scale none`` turns the scaling off
   completely.

4. **Match the background levels** (``--match_background``, optional).  The
   background of an image can be biased by an extended object, which leaves
   steps at the edges of the images in the mosaic.  This option solves for
   additive offsets from the median differences in the overlaps of the images.

5. **Reproject and combine.**  The images are reprojected onto the output WCS
   with ``reproject`` (``--reproject``), and combined with a weighted mean, with
   one weight per image (the inverse of its median variance).  Masked pixels
   never contribute.  Where at least three images contribute, outliers (e.g.,
   cosmic rays not masked by the reduction) are rejected: at each pixel, the
   image that deviates most from the mean of the others is rejected if its
   deviation is larger than ``--rej`` times its error, which includes a
   fractional term (``--rej_frac``, 0.1) so that small differences in the cores
   of bright sources are not rejected.

The errors
----------

The ``UNCERT`` extension of a mosaic is a *noise density*.  The resampling of
the input images leads to correlations between neighboring pixels.  To mitigate
this, the error of a single pixel in the mosaic is overestimated, such that
summing the variance over an aperture much larger than the matching kernel gives
the error of the summed flux (and an upper bound for smaller apertures).  The
errors do not include the errors of the flux scales, which are given in the
``INPUTS`` table.

Execution
---------

An example call looks like this:

.. code-block:: console

   pipnick_mosaic n3_2x2_fast_full.ecsv n3_2x2_slow_full.ecsv --outdir mosaics

The first set of arguments is one or more per-setup metadata files that include
the frames in the mosaic.  At least one file must be provided.

Additional options are:

- ``--outdir``: The output directory.  If not provided, this is set to a
  ``mosaic`` subdirectory of the reduction directory.  If providing more than
  one metadata table, this must always be defined.

- ``--mosaics``: A list of mosaic labels, indicating only those in the list
  should be constructed.

- ``--psf_match``: ``mosaic`` (default), ``group``, or ``none``; see above.

- ``--kernel``: The PSF matching kernel, ``moffat`` (default) or ``gaussian``.

- ``--psf_margin``: The target FWHM in units of the worst FWHM (default 1.05).

- ``--scale``: ``sources`` (default), ``calib``, or ``none``.

- ``--reference``: The reference image of the flux scale of its mosaic.

- ``--match_radius``: The radius (arcsec) for matching sources between images
  (default 1).

- ``--min_common``: The minimum number of sources common to a pair of images
  used for the scaling (default 3).

- ``--min_snr``: The minimum signal-to-noise ratio of the sources used for the
  scaling (default 20).

- ``--background``: ``constant`` (default), ``2d``, or ``none``.

- ``--match_background``: Match the background levels in the overlaps.

- ``--reproject``: ``adaptive`` (default), ``interp``, or ``exact`` (slow).

- ``--pixel_scale``: The output pixel scale (arcsec).

- ``--rej``, ``--rej_frac``, ``--no_reject``: The outlier rejection; see above.

- ``-o``, ``--overwrite``: Rebuild the mosaics, even if they are up to date.

- ``--no_qa``: Do not write the QA figure of each mosaic to ``qa/<mosaic>.png``:
  the mosaic, its coverage, the magnitude offsets of the images, and their PSF
  FWHM before and after the matching.  By default, the figure is redrawn when a
  mosaic is constructed, and drawn for a mosaic that is up to date but lacks
  one.

When rerun, a mosaic is rebuilt only if its images, the settings, or the
membership of its WCS group have changed.  A change in the membership of a WCS
group changes its WCS, so all of its mosaics are rebuilt.

Outputs
-------

The products are written to the output directory (by default, the ``mosaic``
subdirectory of the reduction directory).  Each mosaic is written to
``<mosaic>.fits``, with the following extensions:

- ``PRIMARY``: the mosaic image in electrons per second, on the flux scale of
  the reference image (0 where there are no data);

- ``UNCERT``: its error (a noise density; see above);

- ``MASK``: 1 where there are no data;

- ``EXPMAP``: the total exposure time;

- ``COVERAGE``: the number of images combined for each pixel;

- ``NREJECT``: the number of images whose pixels were rejected;

- ``INPUTS``: a table of the input images, with their tables, setups, nights,
  provenance hashes, exposure times, airmasses, PSF FWHM before and after the
  matching (``fwhm_in``, along the major axis of an elliptical PSF;
  ``fwhm_out``), kernels, the accuracy of the matching (``match_ee``, the ratio
  of the encircled energy of the matched PSF model to that of the target within
  1.5 times the target FWHM), magnitude offsets and scales, the number of
  magnitude differences used, the reference, the relative weights, and the
  background offsets applied by ``--match_background`` (``level``; 0 without
  this option).  Values that are unknown (e.g., an airmass) or
  not measured (e.g., ``fwhm_out`` and ``match_ee`` without PSF matching) are
  masked.  Masked values are stored as NaN, the FITS null value for
  floating-point columns, and :meth:`astropy.table.Table.read` masks them again.

The header records the WCS, the labels (``PNMOSAIC``, ``PNWCSGRP``), the filter
(``mixed`` if the images have different filters), the night, the target PSF, the
settings, the reference image, the number of images on the flux scale
(``PNNSCALE``), and the number of rejected pixels.  Additionally, the header
provides the pipnick version (``PNVERS``) and the hashes that decide whether the
mosaic is up to date (``PNMOSHSH``, ``PNGRPHSH``).  If the reference image has a
photometric calibration, its zero point is given as ``PHOTZP``.  Warnings are
recorded as ``PNWARN`` (images with different filters were combined) and
``PNWARN2`` (images that could not be scaled).

The output directory also has:

- ``<wcs_group>_wcs.fits``: the output WCS of each WCS group (a header only);

- ``summary.ecsv``: a table with one row per mosaic providing the WCS group,
  filter, number of images and of scaled images, target FWHM, status
  (``constructed``, ``up to date``, or ``failed``, with the ``reason``), and
  file;

- ``qa/<mosaic>.png``: the QA figure of each mosaic (see above).

The log is written to the ``log`` subdirectory of the reduction directory
(or of ``--outdir``); see :ref:`common`.

Troubleshooting
---------------

- **"No mosaics are defined"**: The ``mosaic`` columns of the tables are empty
  and must be filled.  Assign them using ``pipnick_metadata raw n3 --split
  --append --assign_mosaics`` or edit the tables directly.

- **"... has no astrometric solution"**, **"... has no source table"**, or
  **"... has no PSF model"**:  Every image assigned to a mosaic must have both
  astrometry and photometry (run ``pipnick_astrometry`` and
  ``pipnick_photometry``).  If one image has no astrometric solution or source
  table, construction of all mosaics stops.  If one image has no PSF model, only
  its mosaic fails (or its whole WCS group, with ``--psf_match group``).  You can
  either run the necessary code or clear the image's ``mosaic`` and
  ``wcs_group`` labels in the table.

- **"Mosaic ... combines images from several nights"**: Mosaics of several
  nights are not yet supported; give the images of each night their own
  labels.

- **"The PSF FWHM of ... is much larger than the median"**: An image in the
  mosaic has particularly poor seeing, which dominates the resolution set for
  the whole mosaic.  You can either remove it from the mosaic (clear its
  ``mosaic`` label) or accept the lower resolution.

- **"The Moffat matching kernel ... failed; using a Gaussian kernel instead"**:
  The PSF of the image is not matched well.  Check ``match_ee`` in the
  ``INPUTS`` table (close to 1 for a good match).

- **Images that are not scaled** (``scaled`` is false in ``INPUTS`` and
  ``PNWARN2``):  The images share fewer than ``--min_common`` sources with the
  other images in the mosaic.  Try lowering ``--min_snr`` or ``--min_common``,
  or use ``--scale calib``.

- **Steps at the edges of the images** in the mosaic:  The backgrounds of the
  images differ, e.g., because an extended object biases them.  Use
  ``--match_background`` and avoid ``--background 2d`` for objects that fill
  much of the field.

- **A failed WCS group** fails all of its mosaics, with the reason in the
  summary; the other groups are still constructed.

Using the pipeline in Python
----------------------------

.. code-block:: python

   from pipnick.pipelines.mosaic import mosaic_all

   files = mosaic_all(
       ['n3_2x2_fast_full.ecsv', 'n3_2x2_slow_full.ecsv'], outdir='mosaics'
   )

See :func:`~pipnick.pipelines.mosaic.mosaic_all` for the full set of
parameters, and :mod:`pipnick.mosaic` for the functions used to match the PSFs,
scale the fluxes, and combine the images.
