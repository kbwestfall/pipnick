"""
QA figure of the sources of an image, e.g., for the astrometry and the
photometry.
"""

import numpy as np

from matplotlib.figure import Figure
from matplotlib.patches import Circle

from pipnick.qa.display import save_figure, show_image


def plot_sources(data, mask, sources, ofile, matched=None, title=None, vector_scale=20.):
    """
    Plot an image with its sources and, if available, the sources matched to
    a reference (e.g., by an astrometric solution) with their residuals.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The image data.
    mask : :class:`numpy.ndarray`
        Boolean array selecting the masked pixels, which are not shown.
    sources : :class:`astropy.table.Table`
        The detected sources, with 0-based pixel coordinates ``x`` and ``y``;
        e.g., see :func:`~pipnick.astrometry.detect.detect_sources`.  They are
        marked by open circles and numbered in order.
    ofile : str, :class:`pathlib.Path`
        The output file.
    matched : dict, optional
        The matched sources: their pixel coordinates (``x`` and ``y``), the
        positions of their matches predicted by the solution (``xref`` and
        ``yref``), and a label for the reference (``label``).  The matched
        sources are marked by crosses, with their residuals shown as vectors.
    title : str, optional
        The figure title.
    vector_scale : float, optional
        The factor used to magnify the residual vectors.
    """
    fig = Figure(figsize=(10, 10))
    ax = fig.add_subplot()
    show_image(ax, data, mask=mask)
    radius = 0.012 * max(data.shape)
    for i, (x, y) in enumerate(zip(sources['x'], sources['y'])):
        ax.add_patch(Circle((x, y), radius, fill=False, color='C0', lw=1))
        ax.text(x + radius, y + radius, str(i + 1), color='C0', fontsize=6)
    if matched is not None and len(matched['x']) > 0:
        ax.scatter(matched['x'], matched['y'], marker='+', color='C1', s=60, lw=1)
        ax.quiver(
            matched['x'], matched['y'],
            vector_scale * (np.asarray(matched['xref']) - np.asarray(matched['x'])),
            vector_scale * (np.asarray(matched['yref']) - np.asarray(matched['y'])),
            angles='xy', scale_units='xy', scale=1, color='C3', width=0.002
        )
        ax.text(
            0.01, 0.99,
            f'{len(matched["x"])} sources (+) matched to {matched["label"]}; residuals (red) '
            f'magnified {vector_scale:g}x',
            transform=ax.transAxes, va='top', color='w', fontsize=9,
            bbox={'facecolor': 'k', 'alpha': 0.5}
        )
    ax.set_xlim(-0.5, data.shape[1] - 0.5)
    ax.set_ylim(-0.5, data.shape[0] - 0.5)
    ax.set_xlabel('x (pixel, 0-based)')
    ax.set_ylabel('y (pixel, 0-based)')
    if title is not None:
        ax.set_title(title)
    fig.tight_layout()
    save_figure(fig, ofile, dpi=100)
