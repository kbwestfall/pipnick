"""
Three-color (RGB) images.

The three images must have the same WCS and shape, i.e., they must be
sampled on the same pixel grid; use ``pipnick_mosaic`` with the images (or
mosaics) in one WCS group to reproject them onto a common grid (see
:mod:`pipnick.pipelines.mosaic`).  The color image is made with
:func:`astropy.visualization.make_lupton_rgb` (the scheme of Lupton et al.
2004) or :func:`astropy.visualization.make_rgb`.
"""

from pathlib import Path

import matplotlib.image
import numpy as np

from astropy.stats import sigma_clipped_stats
from astropy.visualization import LinearStretch, LogStretch, ManualInterval
from astropy.visualization import make_lupton_rgb, make_rgb

from pipnick.astrometry.fit import wcs_from_header
from pipnick.image.io import read_reduced
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickParameterError


def read_image(path):
    """
    Read an image, its good-pixel map, and its WCS.

    Pixels are good if they are not masked (by the ``MASK`` extension, if it
    exists) and finite.

    Parameters
    ----------
    path : str, :class:`pathlib.Path`
        The FITS file, with the image in the primary extension.

    Returns
    -------
    data : :class:`numpy.ndarray`
        The image, 0 at the bad pixels.
    good : :class:`numpy.ndarray`
        Boolean array selecting the good pixels.
    wcs : :class:`astropy.wcs.WCS`
        The celestial WCS, or None if the header has none.
    """
    data, mask, _, header = read_reduced(path)
    good = np.isfinite(data) & np.logical_not(mask)
    return np.where(good, data, 0.), good, wcs_from_header(header)


def check_same_wcs(paths, wcs_list, shapes, tol=0.01):
    """
    Check that images have the same WCS and shape.

    The WCSs are compared by mapping a grid of pixels of the first image to
    the sky and back to the pixels of each of the others, so that equivalent
    WCSs written differently are accepted.

    Parameters
    ----------
    paths : list
        The image files, used in messages.
    wcs_list : list
        The :class:`astropy.wcs.WCS` of each image (None if it has none).
    shapes : list
        The shape of each image.
    tol : float, optional
        The largest acceptable difference in pixels.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the images do not have the same WCS and shape.
    """
    advice = (
        '  Construct the images with pipnick_mosaic, with all of them in the same WCS group '
        '(the wcs_group column of the metadata table), so that they share a WCS; a mosaic can '
        'have a single image.'
    )
    if any(w is None for w in wcs_list):
        missing = [Path(p).name for p, w in zip(paths, wcs_list) if w is None]
        raise PipnickDataModelError(f'{", ".join(missing)} have no WCS.' + advice)
    if len(set(tuple(s) for s in shapes)) > 1:
        raise PipnickDataModelError(
            'The images have different shapes: '
            + ', '.join(f'{Path(p).name} {s[1]}x{s[0]}' for p, s in zip(paths, shapes)) + '.'
            + advice
        )
    ny, nx = shapes[0]
    y, x = np.mgrid[0:ny:max(ny // 10, 1), 0:nx:max(nx // 10, 1)]
    x = np.append(x.ravel(), nx - 1)
    y = np.append(y.ravel(), ny - 1)
    sky = wcs_list[0].pixel_to_world(x, y)
    for path, wcs in zip(paths[1:], wcs_list[1:]):
        xo, yo = wcs.world_to_pixel(sky)
        offset = float(np.max(np.hypot(xo - x, yo - y)))
        if offset > tol:
            raise PipnickDataModelError(
                f'{Path(path).name} does not have the same WCS as {Path(paths[0]).name} (the '
                f'pixels differ by up to {offset:.3g}).' + advice
            )


def normalize_channels(images, goods, percentile=99.5, common=False, subtract=True):
    """
    Normalize the images used for the three channels.

    Each image has its background (the sigma-clipped median of its good
    pixels) subtracted, if ``subtract`` is True, and is divided by the
    ``percentile`` of its good pixels, so that the channels have comparable
    scales, or, if ``common`` is True, all of the images are divided by the
    largest of these percentiles, which keeps their relative fluxes.  The bad
    pixels are 0.

    Returns
    -------
    list
        The normalized images.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if an image has no good pixels, or a percentile is not positive.
    """
    _images = []
    for data, good in zip(images, goods):
        if not np.any(good):
            raise PipnickDataModelError('An image has no good pixels.')
        background = sigma_clipped_stats(data[good], sigma=3.)[1] if subtract else 0.
        _images.append(np.where(good, data - background, 0.))
    levels = [float(np.percentile(d[g], percentile)) for d, g in zip(_images, goods)]
    if any(level <= 0 for level in levels):
        raise PipnickDataModelError(
            f'The {percentile} percentile of an image is not positive; check the images or '
            'the percentile.'
        )
    if common:
        levels = [max(levels)] * 3
    return [d / level for d, level in zip(_images, levels)]


def color_image(
    red, green, blue, method='lupton', stretch=0.5, Q=8., log_a=1000., scales=(1., 1., 1.),
    percentile=99.5, common=False, subtract=True
):
    """
    Construct a three-color image.

    The images are normalized (see :func:`normalize_channels`), multiplied by
    ``scales`` (to adjust the color balance), and combined:

        - ``'lupton'``: with :func:`astropy.visualization.make_lupton_rgb`,
          with ``stretch`` and ``Q``;
        - ``'linear'``: with :func:`astropy.visualization.make_rgb` and a
          linear stretch, from 0 to 1;
        - ``'log'``: with :func:`astropy.visualization.make_rgb` and a
          logarithmic stretch with ``log_a``, from 0 to 1.

    Parameters
    ----------
    red, green, blue : tuple
        The image and good-pixel map of each channel; see :func:`read_image`.
    method : str, optional
        ``'lupton'``, ``'linear'``, or ``'log'``.
    stretch, Q : float, optional
        The parameters of the Lupton et al. (2004) scheme.
    log_a : float, optional
        The parameter of the logarithmic stretch; see
        :class:`astropy.visualization.LogStretch`.
    scales : tuple, optional
        Factors applied to the normalized red, green, and blue images.
    percentile, common, subtract : optional
        See :func:`normalize_channels`.

    Returns
    -------
    :class:`numpy.ndarray`
        The image, with shape (rows, columns, 3) and type ``uint8``.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if the method is unknown.
    """
    images = normalize_channels(
        [red[0], green[0], blue[0]], [red[1], green[1], blue[1]], percentile=percentile,
        common=common, subtract=subtract
    )
    r, g, b = [d * s for d, s in zip(images, scales)]
    if method == 'lupton':
        # NOTE: make_lupton_rgb divides the stretched intensity by the
        # intensity everywhere, and only then replaces the result where the
        # intensity is 0 (e.g., masked pixels); the warning for these
        # divisions is suppressed, and the output has no invalid values
        with np.errstate(invalid='ignore', divide='ignore'):
            return make_lupton_rgb(r, g, b, minimum=0., stretch=stretch, Q=Q)
    if method == 'linear':
        return make_rgb(
            r, g, b, interval=ManualInterval(vmin=0., vmax=1.), stretch=LinearStretch()
        )
    if method == 'log':
        return make_rgb(
            r, g, b, interval=ManualInterval(vmin=0., vmax=1.), stretch=LogStretch(a=log_a)
        )
    raise PipnickParameterError(f'Unknown method {method}; must be lupton, linear, or log.')


def write_color_image(ofile, rgb):
    """
    Write a three-color image (e.g., PNG or JPEG, from the extension of
    ``ofile``), with the first row of the images at the bottom, so that
    north is up for the north-up mosaics of ``pipnick_mosaic``.
    """
    matplotlib.image.imsave(ofile, rgb, origin='lower')
