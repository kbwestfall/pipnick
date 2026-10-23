"""
Aperture photometry at the positions of known sources.
"""

import numpy as np

from photutils.aperture import ApertureStats, CircularAperture


def aperture_analysis(data, x, y, radius, local_bkg=None, error=None, mask=None):
    """
    Measure the flux within circular apertures at the provided positions.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The image data.
    x, y : :class:`numpy.ndarray`
        The 0-based pixel coordinates of the aperture centers.
    radius : float
        The radius of the apertures in pixels.
    local_bkg : :class:`numpy.ndarray`, optional
        The local background per pixel at each position, which is subtracted.
        If None, no background is subtracted.
    error : :class:`numpy.ndarray`, optional
        The 1-sigma error in each pixel, used to compute the flux errors.
    mask : :class:`numpy.ndarray`, optional
        Boolean array selecting the pixels to ignore.

    Returns
    -------
    flux : :class:`numpy.ma.MaskedArray`
        The background-subtracted flux in each aperture.  The flux is masked
        if the aperture includes masked pixels or extends beyond the image.
    flux_err : :class:`numpy.ma.MaskedArray`
        The error in the flux, masked like the flux, and also if ``error`` is
        None.
    """
    _x = np.atleast_1d(np.asarray(x, dtype=float))
    _y = np.atleast_1d(np.asarray(y, dtype=float))
    _mask = np.zeros(data.shape, dtype=bool) if mask is None else mask
    apertures = CircularAperture(np.column_stack([_x, _y]), r=radius)
    stats = ApertureStats(data, apertures, error=error, mask=_mask, local_bkg=local_bkg)

    # Apertures with masked pixels, or that extend beyond the image (which
    # count as masked when the mask is padded)
    nmasked = ApertureStats(_mask.astype(float), apertures).sum
    beyond = (
        (_x - radius < -0.5) | (_x + radius > data.shape[1] - 0.5)
        | (_y - radius < -0.5) | (_y + radius > data.shape[0] - 0.5)
    )
    bad = (np.asarray(nmasked) > 0) | beyond | np.logical_not(np.isfinite(stats.sum))

    flux = np.ma.MaskedArray(np.where(bad, 0., stats.sum), mask=bad)
    if error is None:
        flux_err = np.ma.MaskedArray(np.zeros(len(_x)), mask=np.ones(len(_x), dtype=bool))
    else:
        err_bad = bad | np.logical_not(np.isfinite(stats.sum_err))
        flux_err = np.ma.MaskedArray(np.where(err_bad, 0., stats.sum_err), mask=err_bad)
    return flux, flux_err
