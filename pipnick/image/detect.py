"""
Detection of point sources.

Every step that detects sources uses :func:`find_sources`, with its own
settings, and selects from the result with the selection functions of this
module (:func:`brightest`, :func:`isolated`) or its own:

    - the astrometry detects the brightest sources, excluding those at the
      border (see :func:`~pipnick.pipelines.astrometry.astrometry_all`);
    - the PSF photometry detects the stars that model the PSF (see
      :func:`~pipnick.photometry.starfind.select_psf_stars`), and all of the
      sources, through :class:`SourceFinder`, also in the residual images of
      its later iterations (see
      :func:`~pipnick.photometry.psf_photometry.psf_analysis`).
"""

import numpy as np

from astropy.table import Table
from photutils.detection import DAOStarFinder

from pipnick.image.background import global_background
from pipnick.pkg.exceptions import PipnickParameterError


def find_sources(
    data, mask=None, fwhm=3., thresh=5., threshold=None, background='global', sigma=3.,
    masked='background', min_sep=None, exclude_border=False
):
    """
    Detect the point sources in an image.

    The sources are found by :class:`photutils.detection.DAOStarFinder` in the
    background-subtracted image, above a threshold set by the noise of the
    background.

    Parameters
    ----------
    data : :class:`numpy.ndarray`, :class:`numpy.ma.MaskedArray`
        The image.  The mask of a masked array is combined with ``mask``.
    mask : :class:`numpy.ndarray`, optional
        Boolean array selecting the pixels to ignore (e.g., bad or saturated
        pixels).
    fwhm : float, optional
        The expected FWHM of the sources in pixels.
    thresh : float, optional
        The detection threshold, as a multiple of the noise of the
        background.
    threshold : float, optional
        The detection threshold in the units of the image, which overrides
        ``thresh``.  Required if ``background`` is ``'none'``.
    background : str, optional
        The background: ``'global'``, the sigma-clipped median and standard
        deviation of the unmasked pixels (see
        :func:`~pipnick.image.background.global_background`), which are
        subtracted and set the threshold; or ``'none'``, for an image whose
        background has been subtracted.
    sigma : float, optional
        The clipping threshold of the ``'global'`` background.
    masked : str, optional
        The value of the masked pixels in the image searched for sources,
        which affects the sources near them: ``'background'``, the background
        level; ``'zero'``, 0 before the background is subtracted; or
        ``'keep'``, their values in ``data``.
    min_sep : float, optional
        The minimum separation of the sources, in pixels.  If None, the
        default of :class:`photutils.detection.DAOStarFinder` is used.
    exclude_border : bool, optional
        Exclude the sources within half the size of the detection kernel of
        the border of the image.

    Returns
    -------
    :class:`astropy.table.Table`
        The sources, in the order they are found: their 0-based pixel
        coordinates (``x``, ``y``), ``flux``, ``peak``, ``sharpness``,
        ``roundness1``, and ``roundness2``; see
        :class:`photutils.detection.DAOStarFinder`.  The table metadata
        record the background level (``background``), its noise (``noise``;
        None if not measured), and the threshold (``threshold``).  No
        sources are found if the threshold is not positive.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if ``background`` or ``masked`` is unknown, or the threshold
        is not given for an image without a background.
    """
    if background not in ['global', 'none']:
        raise PipnickParameterError(f'Unknown background {background}; must be global or none.')
    if masked not in ['background', 'zero', 'keep']:
        raise PipnickParameterError(
            f'Unknown treatment of masked pixels {masked}; must be background, zero, or keep.'
        )
    _mask = np.zeros(np.shape(data), dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    if isinstance(data, np.ma.MaskedArray):
        _mask = _mask | np.ma.getmaskarray(data)
        _data = np.asarray(data.data, dtype=float)
    else:
        _data = np.asarray(data, dtype=float)
    if masked == 'zero':
        _data = np.where(_mask, 0., _data)

    level = 0.
    noise = None
    if background == 'global':
        level, noise = global_background(_data, mask=_mask, sigma=sigma)
    if threshold is not None:
        _threshold = threshold
    elif noise is None:
        raise PipnickParameterError(
            'The threshold must be given for an image without a background.'
        )
    else:
        _threshold = thresh * noise
    image = _data - level
    if masked == 'background':
        image = np.where(_mask, 0., image)

    sources = None
    if _threshold > 0:
        finder = DAOStarFinder(
            fwhm=fwhm, threshold=_threshold, exclude_border=exclude_border,
            min_separation=min_sep
        )
        sources = finder(image, mask=_mask)
    columns = ['x', 'y', 'flux', 'peak', 'sharpness', 'roundness1', 'roundness2']
    tbl = Table(
        {name: np.zeros(0) for name in columns} if sources is None else {
            'x': np.asarray(sources['x_centroid'], dtype=float),
            'y': np.asarray(sources['y_centroid'], dtype=float),
            **{n: np.asarray(sources[n], dtype=float) for n in columns[2:]}
        }
    )
    tbl.meta = {'background': float(level), 'noise': noise, 'threshold': float(_threshold)}
    return tbl


class SourceFinder:
    """
    A source finder for :class:`photutils.psf.IterativePSFPhotometry`.

    The finder detects sources with :func:`find_sources` in an image whose
    background has been subtracted, at a fixed threshold, so that the
    residual images of the later iterations of the photometry are searched
    with the threshold set by the noise of the original image.

    Parameters
    ----------
    fwhm : float
        The expected FWHM of the sources in pixels.
    threshold : float
        The detection threshold in the units of the image.
    min_sep : float, optional
        The minimum separation of the sources, in pixels; see
        :func:`find_sources`.
    """
    def __init__(self, fwhm, threshold, min_sep=None):
        self.fwhm = fwhm
        self.threshold = threshold
        self.min_sep = min_sep

    def __call__(self, data, mask=None):
        """
        Detect the sources in an image, as required by
        :class:`photutils.psf.IterativePSFPhotometry`.

        Parameters
        ----------
        data : :class:`numpy.ndarray`
            The background-subtracted image, or a residual image.  Its masked
            pixels are used as they are.
        mask : :class:`numpy.ndarray`, optional
            Boolean array selecting the pixels to ignore.

        Returns
        -------
        :class:`astropy.table.Table`
            The sources; see :func:`find_sources`.  None if no sources are
            found.
        """
        sources = find_sources(
            data, mask=mask, fwhm=self.fwhm, threshold=self.threshold, background='none',
            masked='keep', min_sep=self.min_sep
        )
        return None if len(sources) == 0 else sources


def brightest(sources, n=None):
    """
    Select the brightest sources.

    Parameters
    ----------
    sources : :class:`astropy.table.Table`
        The sources, with their ``flux``; see :func:`find_sources`.
    n : int, optional
        The number of sources to select.  If None, all of them are returned.

    Returns
    -------
    :class:`astropy.table.Table`
        The (at most) ``n`` brightest sources, ordered by decreasing flux.
    """
    return sources[np.argsort(-np.asarray(sources['flux']))][:n]


def isolated(x, y, flux, min_dist):
    """
    Select the sources that are the brightest of their neighbors.

    Parameters
    ----------
    x, y : :class:`numpy.ndarray`
        The pixel coordinates of the sources.
    flux : :class:`numpy.ma.MaskedArray`
        The fluxes of the sources.  Masked fluxes are ignored when the fluxes
        of neighbors are compared.
    min_dist : float
        The distance in pixels within which sources are neighbors.

    Returns
    -------
    :class:`numpy.ndarray`
        Boolean array selecting the sources whose flux is the largest of
        their neighbors (see :func:`numpy.ma.argmax`); a source without
        neighbors is selected.  Callers should also reject sources with a
        masked flux.
    """
    nsources = len(x)
    i, j = np.triu_indices(nsources, k=1)
    dist = np.sqrt((x[j] - x[i])**2 + (y[j] - y[i])**2)
    too_close = np.identity(nsources).astype(bool)
    too_close[i,j] = dist < min_dist
    too_close[j,i] = too_close[i,j]
    _flux = np.tile(flux, (nsources, 1))
    _flux[np.logical_not(too_close)] = np.ma.masked
    return np.ma.argmax(_flux, axis=1) == np.arange(nsources)
