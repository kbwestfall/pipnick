Mosaics
=======

``pipnick_mosaic`` combines reduced science images of the same area of sky
into mosaics: dithered pointings of an extended object, or repeated exposures
of one field (which registers them).  The mosaics of one object in different
filters can share an identical output WCS and shape, so that they can be
combined into colour images.

The images are combined using their astrometric solutions, and the PSF models
and source tables of the photometry, so run ``pipnick_astrometry`` (see
:doc:`astrometry`) and ``pipnick_photometry`` (see :doc:`photometry`) first.
The workflow is:

.. code-block:: console

   pipnick_metadata -d raw --split --assign_mosaics
   pipnick_reduce -t n3_2x2_fast_full.ecsv
   pipnick_astrometry -t n3_2x2_fast_full.ecsv
   pipnick_photometry -t n3_2x2_fast_full.ecsv
   pipnick_mosaic -t n3_2x2_fast_full.ecsv --qa

Defining the mosaics
--------------------

The mosaics are defined by two columns of the per-setup metadata tables,
with string labels:

- ``mosaic``: the frames with the same label are combined into one mosaic;
- ``wcs_group``: the mosaics with the same label share an output WCS and
  shape.  If it is empty, a mosaic is its own WCS group.

An empty label means the frame is not mosaicked.  The labels apply to the
reduced images: all of the frames of a combination group (``comb_id``) must
have the same labels.  A mosaic must belong to one WCS group, and a warning is
issued if it combines images with different filters.

By default, the columns are empty, i.e., no mosaics are defined.  Either edit
the tables, or assign the mosaics automatically with ``pipnick_metadata
--assign_mosaics``:

- The science frames are grouped by pointing: pointings within
  ``--mosaic_radius`` (by default, half of the field of view) of each other,
  directly or through other pointings, form one WCS group.
- Within each pointing group, the frames with the same filter form one mosaic.
- A filter with a single image in a pointing group is not mosaicked.
- The WCS group is labeled with the most common object name of its frames,
  without a trailing filter name, e.g. ``NGC_7331`` for ``NGC 7331 V``, and
  the mosaic with the WCS group and the filter, e.g. ``NGC_7331_V``.  The
  object names are not used to group the frames, because they are often wrong
  (e.g., observers do not always update the filter in the name).
- Only frames without labels are assigned, so edits are kept.  With
  ``--append``, the labels are assigned even if there are no new frames.

A mosaic can have a single image.  To sample several images on the same pixel
grid, e.g., one image of an object in each of three filters for a colour
image, give each image its own ``mosaic`` label and all of them the same
``wcs_group``; each is then only reprojected onto the WCS of the group (it is
not convolved, unless ``--psf_match group`` is used).

The labels are global: the same label in different tables refers to the same
mosaic.  For example, a field observed with two readout modes has two setups,
and so two tables and two reduction directories; give both tables to
``pipnick_mosaic`` (with ``--outdir``) to combine its images.  When the tables
of several setups are built together, ``--assign_mosaics`` gives the same
pointing and filter the same labels in all of them.

.. note::

    For now, a mosaic can only combine images from one night (UT date).  The
    labels, the support for several tables, and the provenance recorded for
    each image are the basis for mosaics of several nights, which will need
    to be enabled explicitly, so that common labels like ``M11_V`` do not
    combine different nights by accident.

How the mosaics are constructed
-------------------------------

The output WCS of each WCS group is a north-up gnomonic projection in the ICRS
that covers all of the images of the group (in all filters), with the finest
pixel scale of the images (or ``--pixel_scale``); see
:func:`reproject.mosaicking.find_optimal_celestial_wcs`.  It is saved to
``<wcs_group>_wcs.fits`` and reused while the group is unchanged, so every
mosaic of the group has the identical WCS and shape.

Each mosaic is then constructed as follows.

1. **Prepare the images.**  Each image is converted to electrons per second,
   and its background is subtracted (``--background``): by default a constant,
   the median of the pixels without sources; or a smooth 2D background
   (256-pixel boxes), which should not be used for objects that fill much of
   the field.

2. **Match the PSFs** (``--psf_match``).  The seeing and focus can change
   during the observations, so by default every image is convolved to the
   worst seeing of its mosaic.  The target is a circular Moffat function with
   ``--psf_margin`` (1.05) times the largest FWHM of the images, and the
   median of their Moffat power indices.  For an elliptical PSF model
   (``pipnick_photometry --fittype ellip``), the FWHM along its major axis is
   used, so that the PSF is broadened to the target in every direction.  Each image is convolved with the
   kernel that matches its PSF model (from the photometry) to the target,
   computed with :func:`photutils.psf_matching.make_kernel` and tapered to a
   circle; if that fails, a Gaussian kernel is used.  For colour images, use
   ``--psf_match group`` to match all of the mosaics of a WCS group to the
   worst seeing of the group.  ``--psf_match none`` turns the matching off.
   A warning is issued if one image has much worse seeing than the others,
   because it degrades the resolution of the whole mosaic.

3. **Scale the fluxes** (``--scale``).  The transparency, airmass, and
   exposure time can differ between the images, so by default each image is
   scaled to the flux scale of a reference image (the image with the lowest
   airmass, or ``--reference``).  The scales are measured from the sources
   common to overlapping images: their fluxes are measured on the matched
   images, in apertures of 1.5 times the target FWHM with a local background,
   and the magnitude offsets of all of the images are fit together, so that
   images that share no sources with the reference are scaled through the
   images they overlap.  Images that share too few sources with the others
   (``--min_common``) are not scaled, with a warning.  ``--scale calib`` uses
   the photometric calibrations of the images instead (see :doc:`calibrate`),
   and ``--scale none`` turns the scaling off.

4. **Match the background levels** (``--match_background``, optional).  The
   background of an image can be biased by an extended object, which leaves
   steps at the edges of the images in the mosaic.  This option solves for
   additive offsets from the median differences in the overlaps of the
   images.

5. **Reproject and combine.**  The images are reprojected onto the output
   WCS with ``reproject`` (``--reproject``), and combined with a weighted
   mean, with one weight per image (the inverse of its median variance).
   Masked pixels never contribute.  Where at least three images contribute,
   outliers (e.g., cosmic rays not masked by the reduction) are rejected: at
   each pixel, the image that deviates most from the mean of the others is
   rejected if its deviation is larger than ``--rej`` times its error, which
   includes a fractional term (``--rej_frac``, 0.1) so that small differences
   in the cores of bright sources are not rejected.

The errors
----------

The ``UNCERT`` extension of a mosaic is a *noise density*.  The convolution
and the reprojection correlate the noise of neighboring pixels, so the error
of a single pixel is overestimated, but summing the variance over an aperture
much larger than the matching kernel gives the error of the summed flux
(and an upper bound for smaller apertures).  The errors do not include the
errors of the flux scales, which are given in the ``INPUTS`` table.

Using the command-line script
-----------------------------

.. code-block:: console

   pipnick_mosaic -t n3_2x2_fast_full.ecsv n3_2x2_slow_full.ecsv --outdir mosaics --qa

The options are:

- ``-t``, ``--raw_table``: One or more per-setup metadata files (required).
- ``--outdir``: The output directory.  If not provided, the ``mosaic``
  subdirectory of the reduction directory, which requires a single table.
- ``--mosaics``: Only construct these mosaics (labels).
- ``--psf_match``: ``mosaic`` (default), ``group``, or ``none``.
- ``--kernel``: The PSF matching kernel, ``moffat`` (default) or
  ``gaussian``.
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
- ``--qa``: Write a QA figure for each mosaic to ``qa/<mosaic>.png``: the
  mosaic, its coverage, the magnitude offsets of the images, and their PSF
  FWHM before and after the matching.

When rerun, a mosaic is rebuilt only if its images, the settings, or the
membership of its WCS group have changed.  A change in the membership of a WCS
group changes its WCS, so all of its mosaics are rebuilt.

Outputs
-------

Each mosaic is written to ``<mosaic>.fits``, with the extensions:

- the primary: the mosaic in electrons per second, on the flux scale of the
  reference image (0 where there are no data);
- ``UNCERT``: its error (a noise density; see above);
- ``MASK``: 1 where there are no data;
- ``EXPMAP``: the total exposure time;
- ``COVERAGE``: the number of images combined;
- ``NREJECT``: the number of images whose pixels were rejected;
- ``INPUTS``: a table of the input images, with their tables, setups, nights,
  provenance hashes, exposure times, airmasses, PSF FWHM before and after the
  matching (``fwhm_in``, along the major axis of an elliptical PSF;
  ``fwhm_out``), kernels, the accuracy of the matching
  (``match_ee``, the ratio of the encircled energy of the matched PSF model
  to that of the target within 1.5 times the target FWHM), magnitude offsets
  and scales, the number of magnitude differences used, the reference, the
  relative weights, and the background levels subtracted.

The header records the labels, the filter, the night, the target PSF, the
settings, the reference image, and the number of rejected pixels.  If the
reference image has a photometric calibration, its zero point is given as
``PHOTZP``.  A summary of the mosaics is written to ``summary.ecsv``.

Colour images
-------------

``pipnick_rgb`` makes a three-colour image from three images with the same
WCS and shape, e.g., the mosaics of one WCS group:

.. code-block:: console

   pipnick_rgb M11_R.fits M11_V.fits M11_B.fits -o M11.png

The images are given in the order red, green, blue.  If they do not have the
same WCS and shape, ``pipnick_rgb`` stops with an error; construct them with
``pipnick_mosaic``, with all of them in the same WCS group.  For colours that
are not affected by differences in seeing, use ``--psf_match group``.

Each image has its background (the median) subtracted (unless
``--no_subtract``) and is divided by its ``--percentile`` (99.5) percentile,
so that the channels have comparable brightness, or, with ``--common``, all of
the images are divided by the largest of these percentiles, which keeps their
relative fluxes.  ``--scales`` multiplies the three normalized images, to
adjust the colour balance.  The image is then made with the scheme of Lupton
et al. (2004; ``--method lupton``, the default, with ``--stretch`` and
``--Q``; see :func:`astropy.visualization.make_lupton_rgb`), or a linear or
logarithmic stretch (``--method linear`` or ``log``, with ``--log_a``; see
:func:`astropy.visualization.make_rgb`); see the `astropy documentation
<https://docs.astropy.org/en/stable/visualization/rgb.html>`__.  Masked
pixels are black.  The format of the output is set by its extension (e.g.,
``png`` or ``jpg``), and north is up for the mosaics of ``pipnick_mosaic``.

Using the pipeline in Python
----------------------------

.. code-block:: python

   from pipnick.pipelines.mosaic import mosaic_all

   files = mosaic_all(['n3_2x2_fast_full.ecsv', 'n3_2x2_slow_full.ecsv'], outdir='mosaics')

See :func:`~pipnick.pipelines.mosaic.mosaic_all` for the full set of
parameters, and :mod:`pipnick.mosaic` for the functions used to match the PSFs,
scale the fluxes, and combine the images.
