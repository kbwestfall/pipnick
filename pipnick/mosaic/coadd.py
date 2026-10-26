r"""
Reprojection and combination of images into a mosaic.

The images are reprojected onto the output WCS with ``reproject``.  So that
masked pixels never contribute (and are never represented by NaN), each image
:math:`d` is reprojected as the product :math:`d\,g` with its good-pixel map
:math:`g`, together with :math:`g` itself and its variance times :math:`g`;
the reprojected image is :math:`(d\,g)' / g'`, and pixels with :math:`g'`
below a threshold are masked.  Pixels outside the footprint of an image, which
``reproject`` returns as NaN, are masked immediately.

The values are per pixel, so they are multiplied by the ratio of the output
and input pixel areas.  The variance is multiplied by the same ratio: like the
convolution (see :mod:`pipnick.mosaic.psfmatch`), the reprojection correlates
the noise of neighboring pixels, and the reprojected variance is a variance
density, correct for sums over apertures much larger than a pixel.

The reprojected images are combined with one weight per image, the inverse
of its typical (median) variance, rather than per-pixel inverse-variance
weights, which include the Poisson noise of the sources and would bias bright
sources low; see :func:`combine`.
"""

import numpy as np

from astropy import units
from astropy.stats import sigma_clipped_stats
from astropy.wcs.utils import proj_plane_pixel_area
from reproject import reproject_adaptive, reproject_exact, reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs

from pipnick.pkg.exceptions import PipnickParameterError


def output_wcs(wcs_list, shape_list, pixel_scale=None):
    """
    Find the output WCS and shape that cover a set of images.

    The WCS is a north-up gnomonic (TAN) projection in the ICRS; see
    :func:`reproject.mosaicking.find_optimal_celestial_wcs`.

    Parameters
    ----------
    wcs_list : list
        The :class:`astropy.wcs.WCS` of each image.
    shape_list : list
        The shape of each image.
    pixel_scale : float, optional
        The output pixel scale (arcsec).  If None, the finest pixel scale of
        the images is used.

    Returns
    -------
    wcs : :class:`astropy.wcs.WCS`
        The output WCS.
    shape : tuple
        The output shape.
    """
    if pixel_scale is None:
        _scale = min(np.sqrt(proj_plane_pixel_area(w)) * 3600 for w in wcs_list)
    else:
        _scale = float(pixel_scale)
    wcs, shape = find_optimal_celestial_wcs(
        [(shape, w) for shape, w in zip(shape_list, wcs_list)], frame='icrs',
        resolution=_scale * units.arcsec
    )
    return wcs, tuple(int(s) for s in shape)


def footprint_bbox(wcs_in, shape_in, wcs_out, shape_out, pad=2):
    """
    Return the bounding box of an image on the output pixel grid.

    Parameters
    ----------
    wcs_in : :class:`astropy.wcs.WCS`
        The WCS of the image.
    shape_in : tuple
        The shape of the image.
    wcs_out : :class:`astropy.wcs.WCS`
        The output WCS.
    shape_out : tuple
        The output shape.
    pad : int, optional
        Padding added on each side (pixels).

    Returns
    -------
    tuple
        The slices ``(rows, columns)`` of the bounding box, or None if the
        image does not overlap the output grid.
    """
    ny, nx = shape_in
    # The edges of the image (pixel corners), sampled finely enough to
    # capture distortion
    t = np.linspace(0., 1., 21)
    x = -0.5 + t * nx
    y = -0.5 + t * ny
    xs = np.concatenate([x, np.full(21, nx - 0.5), x, np.full(21, -0.5)])
    ys = np.concatenate([np.full(21, -0.5), y, np.full(21, ny - 0.5), y])
    xo, yo = wcs_out.world_to_pixel(wcs_in.pixel_to_world(xs, ys))
    x0 = max(int(np.floor(np.min(xo))) - pad, 0)
    x1 = min(int(np.ceil(np.max(xo))) + pad + 1, shape_out[1])
    y0 = max(int(np.floor(np.min(yo))) - pad, 0)
    y1 = min(int(np.ceil(np.max(yo))) + pad + 1, shape_out[0])
    if x1 <= x0 or y1 <= y0:
        return None
    return slice(y0, y1), slice(x0, x1)


def reproject_image(
    data, var, good, wcs_in, wcs_out, shape_out, method='adaptive', min_good=0.99
):
    r"""
    Reproject an image, its variance, and its good-pixel map.

    Parameters
    ----------
    data, var : :class:`numpy.ndarray`
        The image and its variance (0 at masked pixels).
    good : :class:`numpy.ndarray`
        Boolean array selecting the good pixels.
    wcs_in : :class:`astropy.wcs.WCS`
        The WCS of the image.
    wcs_out : :class:`astropy.wcs.WCS`
        The output WCS (e.g., of a bounding box; see :func:`footprint_bbox`).
    shape_out : tuple
        The output shape.
    method : str, optional
        The reprojection function: ``'adaptive'``
        (:func:`reproject.reproject_adaptive`), ``'interp'``
        (:func:`reproject.reproject_interp`, bilinear), or ``'exact'``
        (:func:`reproject.reproject_exact`, slow).
    min_good : float, optional
        Output pixels with a reprojected good-pixel fraction below this are
        masked.

    Returns
    -------
    data, var : :class:`numpy.ndarray`
        The reprojected image and variance, per output pixel, 0 at the masked
        pixels.
    good : :class:`numpy.ndarray`
        Boolean array selecting the good output pixels.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if the method is unknown.
    """
    if method == 'adaptive':
        function = reproject_adaptive
        kwargs = {'conserve_flux': False, 'boundary_mode': 'ignore'}
    elif method == 'interp':
        function = reproject_interp
        kwargs = {'order': 'bilinear'}
    elif method == 'exact':
        function = reproject_exact
        kwargs = {}
    else:
        raise PipnickParameterError(f'Unknown reprojection method {method}.')
    g = np.asarray(good, dtype=float)
    planes = []
    for plane in [np.where(good, data, 0.), np.where(good, var, 0.), g]:
        array, footprint = function((plane, wcs_in), wcs_out, shape_out=shape_out, **kwargs)
        # NaN marks pixels outside the image; they are masked immediately
        covered = (footprint > 0) & np.isfinite(array)
        planes.append((np.where(covered, array, 0.), covered))
    (dg, covered), (vg, _), (gr, _) = planes
    _good = covered & (gr >= min_good)
    _gr = np.where(_good, gr, 1.)
    ratio = proj_plane_pixel_area(wcs_out) / proj_plane_pixel_area(wcs_in)
    _data = np.where(_good, dg / _gr, 0.) * ratio
    _var = np.where(_good, np.maximum(vg / _gr, 0.), 0.) * ratio
    return _data, _var, _good


def combine(images, shape, weights, rej=3., rej_frac=0.1, maxiter=3, min_reject=3):
    r"""
    Combine reprojected images.

    The combined image is the weighted mean

    .. math::

        M = \frac{\sum_i w_i g_i d_i}{\sum_i w_i g_i},\qquad
        V = \frac{\sum_i w_i^2 g_i v_i}{(\sum_i w_i g_i)^2},

    where :math:`w_i` is the weight of image :math:`i`, :math:`g_i` selects
    its used pixels, and :math:`v_i` is its variance.

    Where at least ``min_reject`` images contribute, outliers (e.g., cosmic
    rays and other artifacts not masked by the reduction) are rejected.  At
    each pixel, the image that deviates most from the mean of the *other*
    images, :math:`M_{-i}`, is rejected if its deviation is larger than
    ``rej`` times :math:`\sqrt{v_i + V_{-i} + (f M_{-i})^2}`, where :math:`f`
    is ``rej_frac``.  The fractional term allows for small differences between
    the images that are significant compared to the photon noise of bright
    sources, e.g. from the PSF matching, flux scaling, and registration, so
    that the cores of bright sources are not rejected.  Comparing to the mean of the
    other images (and rejecting only one image per pixel and iteration)
    prevents an outlier from causing the rejection of the good pixels of the
    other images.  This is repeated up to ``maxiter`` times.

    Parameters
    ----------
    images : list
        For each image, a tuple ``(bbox, data, var, good)``, with the slices
        of its bounding box on the output grid (see :func:`footprint_bbox`)
        and its reprojected image, variance, and good-pixel map (see
        :func:`reproject_image`) within the bounding box.
    shape : tuple
        The output shape.
    weights : array-like
        The weight of each image.
    rej : float, optional
        The rejection threshold.  If None, no pixels are rejected.
    rej_frac : float, optional
        The fractional error added to the rejection threshold; see above.
    maxiter : int, optional
        The maximum number of rejection iterations.
    min_reject : int, optional
        The minimum number of contributing images for rejection.

    Returns
    -------
    data, var : :class:`numpy.ndarray`
        The combined image and its variance, 0 where no image contributes.
    coverage : :class:`numpy.ndarray`
        The number of images contributing to each pixel.
    used : list
        For each image, the boolean array of its pixels used in the
        combination (its good pixels that were not rejected).
    nrej : int
        The number of rejected pixels.
    """
    used = [np.asarray(g, dtype=bool).copy() for _, _, _, g in images]
    nrej = 0
    for iteration in range(maxiter + 1):
        sw, swd, sw2v, coverage = _sums(images, used, weights, shape)
        if rej is None or iteration == maxiter:
            break
        # The largest normalized deviation from the mean of the other images,
        # and the image with it, at each pixel
        worst = np.zeros(shape, dtype=float)
        which = np.full(shape, -1, dtype=int)
        for k, ((bbox, d, v, _), u, w) in enumerate(zip(images, used, weights)):
            rest = sw[bbox] - w * u
            ok = u & (rest > 0) & (coverage[bbox] >= min_reject)
            _rest = np.where(ok, rest, 1.)
            mean = (swd[bbox] - w * u * d) / _rest
            var = (sw2v[bbox] - w**2 * u * v) / _rest**2
            err = np.sqrt(np.maximum(v + var, 0.) + (rej_frac * mean)**2)
            dev = np.where(ok & (err > 0), np.abs(d - mean) / np.where(err > 0, err, 1.), 0.)
            update = dev > worst[bbox]
            worst[bbox] = np.where(update, dev, worst[bbox])
            which[bbox] = np.where(update, k, which[bbox])
        new = 0
        for k, ((bbox, _, _, _), u) in enumerate(zip(images, used)):
            deviant = u & (which[bbox] == k) & (worst[bbox] > rej)
            u &= np.logical_not(deviant)
            new += int(np.sum(deviant))
        nrej += new
        if new == 0:
            break
    covered = sw > 0
    _sw = np.where(covered, sw, 1.)
    mean = np.where(covered, swd / _sw, 0.)
    var = np.where(covered, sw2v / _sw**2, 0.)
    return mean, var, coverage, used, nrej


def _sums(images, used, weights, shape):
    """
    Return the sums of the weights, the weighted data, and the squared
    weights times the variances of the used pixels of the images, and the
    number of images used at each pixel; see :func:`combine`.
    """
    sw = np.zeros(shape, dtype=float)
    swd = np.zeros(shape, dtype=float)
    sw2v = np.zeros(shape, dtype=float)
    coverage = np.zeros(shape, dtype=int)
    for (bbox, d, v, _), u, w in zip(images, used, weights):
        sw[bbox] += w * u
        swd[bbox] += w * u * d
        sw2v[bbox] += w**2 * u * v
        coverage[bbox] += u
    return sw, swd, sw2v, coverage


def overlap_differences(images, min_pixels=1000):
    """
    Measure the differences between the levels of overlapping reprojected
    images.

    For each pair of images that overlap in at least ``min_pixels`` good
    pixels, the difference is the sigma-clipped median of the pixel
    differences in the overlap, which rejects the sources.  Its error is the
    clipped standard deviation divided by the square root of the number of
    pixels; this ignores the correlation of the noise of neighboring pixels,
    but only sets the relative weights of the pairs.

    Parameters
    ----------
    images : list
        For each image, a tuple ``(bbox, data, var, good)``; see
        :func:`combine`.
    min_pixels : int, optional
        The minimum number of overlapping good pixels.

    Returns
    -------
    list
        For each overlapping pair, a tuple ``(i, j, diff, err)`` with the
        difference (image ``i`` minus image ``j``) and its error, as
        one-element arrays; see :func:`~pipnick.mosaic.scaling.solve_offsets`.
    """
    measurements = []
    for i, (bi, di, _, gi) in enumerate(images):
        for j in range(i + 1, len(images)):
            bj, dj, _, gj = images[j]
            rows = slice(max(bi[0].start, bj[0].start), min(bi[0].stop, bj[0].stop))
            cols = slice(max(bi[1].start, bj[1].start), min(bi[1].stop, bj[1].stop))
            if rows.stop <= rows.start or cols.stop <= cols.start:
                continue
            si = (
                slice(rows.start - bi[0].start, rows.stop - bi[0].start),
                slice(cols.start - bi[1].start, cols.stop - bi[1].start)
            )
            sj = (
                slice(rows.start - bj[0].start, rows.stop - bj[0].start),
                slice(cols.start - bj[1].start, cols.stop - bj[1].start)
            )
            both = gi[si] & gj[sj]
            n = int(np.sum(both))
            if n < min_pixels:
                continue
            _, median, std = sigma_clipped_stats((di[si] - dj[sj])[both], sigma=3.)
            measurements.append(
                (i, j, np.array([median]), np.array([max(std, 1e-30) / np.sqrt(n)]))
            )
    return measurements
