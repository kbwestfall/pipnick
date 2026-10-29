"""
The background of an image.

Three estimates of the background are used: a global sigma-clipped level
and noise (:func:`global_background`), e.g., for the detection of sources (see
:func:`~pipnick.image.detect.find_sources`); a smooth, 2D background
(:func:`subtract_background`), e.g., for the PSF photometry; and a constant
level measured with the sources masked by image segmentation
(:func:`masked_median`), e.g., for the mosaics.  :func:`estimate_background`
selects between them.  The local backgrounds measured around sources are
part of the photometry.
"""

import numpy as np

from astropy.stats import sigma_clipped_stats
from photutils.background import Background2D
from photutils.segmentation import detect_sources, detect_threshold

from pipnick.pkg.exceptions import PipnickParameterError


def global_background(data, mask=None, sigma=3.):
    """
    Return the global background level and noise of an image: the
    sigma-clipped median and standard deviation of its unmasked pixels; see
    :func:`astropy.stats.sigma_clipped_stats`.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The image.
    mask : :class:`numpy.ndarray`, optional
        Boolean array selecting the pixels to ignore.
    sigma : float, optional
        The clipping threshold, in standard deviations.

    Returns
    -------
    level, noise : float
        The median and standard deviation of the background.
    """
    _, level, noise = sigma_clipped_stats(data, mask=mask, sigma=sigma)
    return level, noise


def subtract_background(data, mask=None, box_size=64):
    """
    Subtract a smooth, 2D background from an image.

    The background is measured with :class:`photutils.background.Background2D`
    in boxes of ``box_size`` pixels, ignoring masked pixels.  The source
    detection used by the PSF photometry assumes background-subtracted data.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The image data.
    mask : :class:`numpy.ndarray`, optional
        Boolean array selecting the pixels to ignore.
    box_size : int, optional
        The size of the boxes used to measure the background.

    Returns
    -------
    :class:`numpy.ndarray`
        The background-subtracted data.
    """
    bkg = Background2D(data, box_size, mask=mask)
    return data - bkg.background


def masked_median(data, mask=None, nsigma=2., npixels=10, dilate=11):
    """
    Return the background level of an image, measured with its sources
    masked.

    The sources are found by image segmentation (see
    :func:`photutils.segmentation.detect_sources`): groups of at least
    ``npixels`` connected pixels above ``nsigma`` times the noise, dilated by
    a ``dilate``-pixel footprint.  The level is the sigma-clipped median of
    the remaining pixels.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The image.
    mask : :class:`numpy.ndarray`
        Boolean array selecting the pixels to ignore.
    nsigma : float, optional
        The detection threshold of the sources, in units of the noise.
    npixels : int, optional
        The minimum number of connected pixels of a source.
    dilate : int, optional
        The size of the footprint used to dilate the source mask.

    Returns
    -------
    float
        The background level.
    """
    threshold = detect_threshold(data, nsigma, mask=mask)
    segm = detect_sources(data, threshold, npixels, mask=mask)
    sources = (
        np.zeros(data.shape, dtype=bool) if segm is None else segm.make_source_mask(size=dilate)
    )
    _mask = sources if mask is None else mask | sources
    _, median, _ = sigma_clipped_stats(data, mask=_mask, sigma=3.)
    return median


def estimate_background(data, mask=None, method='constant', box_size=64):
    """
    Return the background of an image.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The image.
    mask : :class:`numpy.ndarray`, optional
        Boolean array selecting the pixels to ignore.
    method : str, optional
        ``'constant'``, the median with the sources masked (see
        :func:`masked_median`); ``'2d'``, a smooth background in boxes of
        ``box_size`` pixels (see :func:`subtract_background`); or ``'none'``,
        no background (0).
    box_size : int, optional
        The size of the boxes of the ``'2d'`` background.

    Returns
    -------
    float, :class:`numpy.ndarray`
        The background: a constant, or an array with the shape of the image
        for ``'2d'``.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if the method is unknown.
    """
    if method == 'none':
        return 0.
    if method == 'constant':
        return masked_median(data, mask=mask)
    if method == '2d':
        return Background2D(data, box_size, mask=mask).background
    raise PipnickParameterError(
        f'Unknown background method {method}; options are constant, 2d, and none.'
    )
