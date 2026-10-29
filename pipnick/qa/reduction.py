"""
QA figure of a reduction product: a combined calibration frame or a reduced
science image.
"""

import numpy as np

from matplotlib.figure import Figure

from pipnick.qa.display import save_figure, show_image


def plot_product(data, mask, ofile, title=None, unit=None):
    """
    Plot a reduction product: the image, and the histogram of its unmasked
    pixel values.

    The histogram is limited to the 0.1 to 99.9 percentiles of the values,
    so that a few extreme pixels (e.g., saturated stars) do not hide the
    distribution of the rest, and its counts are shown on a logarithmic
    scale.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The image.
    mask : :class:`numpy.ndarray`
        Boolean array selecting the masked pixels, which are not shown.
    ofile : str, :class:`pathlib.Path`
        The output file.
    title : str, optional
        The figure title.
    unit : str, optional
        The unit of the pixel values, used to label the histogram.
    """
    fig = Figure(figsize=(13, 5.5))
    ax_img, ax_hist = fig.subplots(1, 2, width_ratios=[1.3, 1])
    im = show_image(ax_img, data, mask=mask)
    fig.colorbar(im, ax=ax_img, fraction=0.046, pad=0.02)
    ax_img.set_xlabel('x (pixel, 0-based)')
    ax_img.set_ylabel('y (pixel, 0-based)')
    good = np.logical_not(np.asarray(mask, dtype=bool))
    values = data[good]
    if values.size > 0:
        lo, hi = np.percentile(values, [0.1, 99.9])
        if hi <= lo:
            hi = lo + 1.
        ax_hist.hist(values, bins=200, range=(lo, hi), color='C0', histtype='stepfilled')
        ax_hist.set_yscale('log')
        ax_hist.axvline(np.median(values), color='k', lw=0.8, label='median')
        ax_hist.legend(loc='upper right')
    ax_hist.set_xlabel('pixel value' if unit is None else f'pixel value ({unit})')
    ax_hist.set_ylabel('number of pixels')
    ax_hist.set_title(f'{np.sum(np.logical_not(good))} masked pixels', fontsize=9)
    if title is not None:
        fig.suptitle(title)
    fig.tight_layout()
    save_figure(fig, ofile, dpi=100)
