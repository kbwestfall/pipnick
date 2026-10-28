"""
QA figure of a mosaic.
"""

import numpy as np

from matplotlib.figure import Figure

from pipnick.qa.display import save_figure, show_image


def plot_mosaic(hdus, inputs, ofile):
    """
    Write the QA figure of a mosaic: the mosaic, its coverage, the magnitude
    offsets of the images, and their PSF FWHM before and after the matching.

    Parameters
    ----------
    hdus : :class:`astropy.io.fits.HDUList`
        The mosaic; see :func:`~pipnick.pipelines.mosaic.mosaic_all`.
    inputs : :class:`astropy.table.Table`
        The table of the input images.
    ofile : str, :class:`pathlib.Path`
        The output file.
    """
    data = hdus[0].data
    good = hdus['MASK'].data == 0
    fig = Figure(figsize=(12, 10))
    axes = fig.subplots(2, 2)
    show_image(axes[0,0], data, mask=np.logical_not(good))
    axes[0,0].set_title(f'{hdus[0].header["PNMOSAIC"]} ({hdus[0].header["FILTER"]})')
    im = axes[0,1].imshow(hdus['COVERAGE'].data, origin='lower', cmap='viridis')
    fig.colorbar(im, ax=axes[0,1], label='number of images')
    axes[0,1].set_title('coverage')
    index = np.arange(len(inputs))
    axes[1,0].errorbar(index, inputs['offset'], yerr=inputs['offset_err'], fmt='o')
    unscaled = np.logical_not(np.asarray(inputs['scaled']))
    axes[1,0].scatter(
        index[unscaled], np.asarray(inputs['offset'])[unscaled], marker='x', color='C3', s=60,
        label='not scaled'
    )
    axes[1,0].set_xlabel('image')
    axes[1,0].set_ylabel('magnitude offset (mag)')
    axes[1,0].set_title(f'reference: {hdus[0].header["PNREFIMG"]}', fontsize=9)
    axes[1,1].plot(index, inputs['fwhm_in'], 'o', label='input')
    fwhm_out = np.asarray(inputs['fwhm_out'])
    axes[1,1].plot(index[fwhm_out > 0], fwhm_out[fwhm_out > 0], 's', label='matched')
    if 'PNTGFWHM' in hdus[0].header:
        axes[1,1].axhline(hdus[0].header['PNTGFWHM'], color='k', lw=0.5, label='target')
    axes[1,1].set_xlabel('image')
    axes[1,1].set_ylabel('PSF FWHM (arcsec)')
    axes[1,1].legend()
    fig.tight_layout()
    save_figure(fig, ofile)
