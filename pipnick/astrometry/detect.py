"""
Detect point sources for astrometric calibration.
"""

import numpy as np

from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from photutils.detection import DAOStarFinder


def detect_sources(data, mask=None, fwhm=3., thresh=5., nmax=100):
    """
    Detect the brightest point sources in an image.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The image data.
    mask : :class:`numpy.ndarray`, optional
        Boolean array selecting pixels to ignore (e.g., bad or saturated
        pixels).
    fwhm : float, optional
        The expected FWHM of the sources in pixels.
    thresh : float, optional
        The detection threshold, as a multiple of the standard deviation of the
        background.
    nmax : int, optional
        The maximum number of sources to return.

    Returns
    -------
    :class:`astropy.table.Table`
        The 0-based pixel coordinates (``x`` and ``y``) and the flux (``flux``)
        of the detected sources, ordered by decreasing flux.
    """
    _mask = np.zeros(data.shape, dtype=bool) if mask is None else mask
    _, median, std = sigma_clipped_stats(data, mask=_mask, sigma=3.)
    sources = None
    if std > 0:
        finder = DAOStarFinder(fwhm=fwhm, threshold=thresh * std, exclude_border=True)
        sources = finder(data - median, mask=_mask)
    if sources is None or len(sources) == 0:
        return Table({'x': np.zeros(0), 'y': np.zeros(0), 'flux': np.zeros(0)})
    sources = sources[np.argsort(-sources['flux'])][:nmax]
    return Table(
        {
            'x': np.asarray(sources['x_centroid'], dtype=float),
            'y': np.asarray(sources['y_centroid'], dtype=float),
            'flux': np.asarray(sources['flux'], dtype=float)
        }
    )
