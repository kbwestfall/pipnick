.. _rgb:

Three-color Image Construction
==============================

``pipnick_rgb`` makes a three-colour image from three images of the same
area of sky, e.g., the mosaics of one object in three filters.

The algorithm
-------------

1. **Check the images.**  The three images must have the same WCS and shape,
   so that their pixels are aligned.
2. **Normalize each image.**  Its background (the median) is subtracted
   (unless ``--no_subtract``), and it is divided by its ``--percentile``
   (99.5) percentile, so that the channels have comparable brightness.  With
   ``--common``, all of the images are instead divided by the largest of these
   percentiles, which keeps their relative fluxes.  ``--scales`` multiplies the
   three normalized images, to adjust the colour balance.
3. **Stretch.**  The image is made with the scheme of Lupton et al. (2004;
   ``--method lupton``, the default, with ``--stretch`` and ``--Q``; see
   :func:`astropy.visualization.make_lupton_rgb`), or a linear or logarithmic
   stretch (``--method linear`` or ``log``, with ``--log_a``; see
   :func:`astropy.visualization.make_rgb`).  See the `astropy documentation
   <https://docs.astropy.org/en/stable/visualization/rgb.html>`__ for the
   effect of each parameter.  Masked pixels are black.

Inputs
------

- **Three images with the same WCS and shape**, given in the order red,
  green, blue; e.g., the mosaics of one WCS group made by ``pipnick_mosaic``
  (see :doc:`mosaic`).  The masks of the images (the ``MASK`` extension) are
  used if they exist.

Running it
----------

.. code-block:: console

    pipnick_rgb M11_R.fits M11_V.fits M11_B.fits -o M11.png

The options are (the logging options are common to every script; see
:doc:`common`):

- ``red``, ``green``, ``blue``: The images for the three channels.
- ``-o``, ``--ofile``: The output image (default ``rgb.png``); its format
  (e.g., ``png`` or ``jpg``) is set by its extension.
- ``--method``: ``lupton`` (default), ``linear``, or ``log``.
- ``--stretch``, ``--Q``: The linear stretch (default 0.5, in units of the
  normalized images) and the softening (default 8) of the Lupton scheme.
- ``--log_a``: The parameter of the logarithmic stretch (default 1000).
- ``--percentile``: The percentile used to normalize each image (default
  99.5).
- ``--common``: Normalize all of the images by the largest of their
  percentiles.
- ``--no_subtract``: Do not subtract the background of each image.
- ``--scales``: Three factors applied to the normalized red, green, and blue
  images (default 1 1 1).

Outputs
-------

- **The colour image** (``--ofile``), with north up for the mosaics of
  ``pipnick_mosaic``.
- **The log**, in the ``log`` subdirectory of the directory of the output
  image; see :doc:`common`.

Troubleshooting
---------------

- **The images do not have the same WCS and shape**: ``pipnick_rgb`` stops
  with an error.  Construct them with ``pipnick_mosaic``, with all of them in
  the same WCS group (a mosaic can have a single image; see :doc:`mosaic`).
- **Coloured halos around stars**: the images have different seeing.  Match
  their PSFs to the worst seeing of the group with ``pipnick_mosaic
  --psf_match group``.
- **Black cores with coloured rims in bright stars**: their saturated
  pixels are masked, and masked pixels are black.  The rims are where the
  masks of the three images differ.
- **One colour dominates**: adjust ``--scales``, or use ``--common`` to keep
  the relative fluxes of the images (e.g., if they are calibrated).

Using the function in Python
----------------------------

.. code-block:: python

    from pipnick.utils.rgb import check_same_wcs, color_image, read_image, write_color_image

    paths = ['M11_R.fits', 'M11_V.fits', 'M11_B.fits']
    images = [read_image(p) for p in paths]
    check_same_wcs(paths, [i[2] for i in images], [i[0].shape for i in images])
    rgb = color_image(*[(data, good) for data, good, _ in images], method='lupton')
    write_color_image('M11.png', rgb)

See :mod:`pipnick.utils.rgb`.
