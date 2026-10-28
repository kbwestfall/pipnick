"""
Shared tools for the QA figures.
"""

from pathlib import Path

import numpy as np

from astropy.visualization import ZScaleInterval


def show_image(ax, data, mask=None, cmap='gray'):
    """
    Show an image, scaled by the ZScale limits of its unmasked pixels.

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        The axes.
    data : :class:`numpy.ndarray`
        The image.
    mask : :class:`numpy.ndarray`, optional
        Boolean array selecting the masked pixels, which are not shown.
    cmap : str, optional
        The colormap.

    Returns
    -------
    :class:`matplotlib.image.AxesImage`
        The displayed image.
    """
    _mask = np.zeros(data.shape, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    good = np.logical_not(_mask)
    vmin, vmax = ZScaleInterval().get_limits(data[good]) if np.any(good) else (0., 1.)
    return ax.imshow(
        np.ma.MaskedArray(data, mask=_mask), origin='lower', cmap=cmap, vmin=vmin, vmax=vmax
    )


def save_figure(fig, ofile, dpi=None):
    """
    Save a figure, creating its directory if necessary.

    Parameters
    ----------
    fig : :class:`matplotlib.figure.Figure`
        The figure.
    ofile : str, :class:`pathlib.Path`
        The output file; the format is set by its extension (e.g., ``.png``
        or ``.pdf``).
    dpi : float, optional
        The resolution; the matplotlib default if None.
    """
    _ofile = Path(ofile)
    _ofile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_ofile, dpi=dpi)
