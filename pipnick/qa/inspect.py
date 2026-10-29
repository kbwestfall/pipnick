"""
Display of an image with its sources, for ``pipnick_inspect``.
"""

import numpy as np

from matplotlib import pyplot
from matplotlib.figure import Figure
from matplotlib.patches import Circle

from pipnick.qa.display import save_figure, show_image


def plot_inspection(data, mask, sources, removed=None, highlight=None, title=None, ofile=None):
    """
    Show an image with its measured sources, the sources removed by the
    photometry, and a highlighted source.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The image.
    mask : :class:`numpy.ndarray`
        Boolean array selecting the masked pixels, which are not shown.
    sources : :class:`astropy.table.Table`
        The measured sources, with their ``id`` and 0-based pixel coordinates
        ``x`` and ``y``; e.g., the source table of the image (see
        :func:`~pipnick.pipelines.photometry.photometry_all`).  They are
        marked by circles, labeled by their ``id``.
    removed : :class:`astropy.table.Table`, optional
        The sources removed by the photometry, with the pixel coordinates
        where they were detected (``x_init``, ``y_init``); see
        :func:`~pipnick.pipelines.photometry.vet_sources`.  They are marked by
        red crosses.
    highlight : tuple, optional
        The pixel coordinates of a source to highlight, e.g., a traced source.
    title : str, optional
        The figure title.
    ofile : str, :class:`pathlib.Path`, optional
        The output file.  If None, the figure is shown in an interactive
        window, which requires a display.
    """
    fig = pyplot.figure(figsize=(10, 10)) if ofile is None else Figure(figsize=(10, 10))
    ax = fig.add_subplot()
    show_image(ax, data, mask=mask)
    radius = 0.012 * max(data.shape)
    for i, x, y in zip(sources['id'], sources['x'], sources['y']):
        ax.add_patch(Circle((x, y), radius, fill=False, color='C0', lw=1))
        ax.text(x + radius, y + radius, str(i), color='C0', fontsize=6)
    if removed is not None and len(removed) > 0:
        ax.scatter(
            np.asarray(removed['x_init']), np.asarray(removed['y_init']), marker='x',
            color='C3', s=30, lw=1, label=f'{len(removed)} removed sources'
        )
        ax.legend(loc='upper right')
    if highlight is not None:
        ax.add_patch(Circle(highlight, 2.5 * radius, fill=False, color='yellow', lw=2))
    ax.set_xlim(-0.5, data.shape[1] - 0.5)
    ax.set_ylim(-0.5, data.shape[0] - 0.5)
    ax.set_xlabel('x (pixel, 0-based)')
    ax.set_ylabel('y (pixel, 0-based)')
    if title is not None:
        ax.set_title(title)
    fig.tight_layout()
    if ofile is None:
        pyplot.show()
        return
    save_figure(fig, ofile, dpi=100)
