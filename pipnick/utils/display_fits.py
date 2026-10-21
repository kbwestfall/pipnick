from pathlib import Path
from typing import Union

import numpy as np
from matplotlib import pyplot, cm

from astropy.io import fits
from astropy.visualization import ZScaleInterval
from astropy.nddata import CCDData

from pipnick.utils.fits_class import Fits_Simple
from pipnick.utils.dir_nav import unzip_directories


def print_fits_info(image_path: str):
    """
    Print HDU List information and display the FITS image data.

    Parameters
    ----------
    image_path : str
        Path to the FITS image (greyscale only).
    """
    with fits.open(image_path) as hdul:
        print("\nHDU Header")
        print(repr(hdul[0].header))
        
        pyplot.figure(figsize=(8, 6))
        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(hdul[0].data)
        pyplot.imshow(hdul[0].data, origin='lower', vmin=vmin, vmax=vmax)
        pyplot.gcf().set_dpi(300)
        pyplot.colorbar()
        pyplot.show()


def display_nickel(image: Union[str, Path, Fits_Simple]):
    """
    Display the data of a FITS image using zscale coloring, with masked
    pixels colored red.

    Parameters
    ----------
    image : str, :class:`pathlib.Path`, :class:`~pipnick.utils.fits_class.Fits_Simple`
        The Fits_Simple object or path to the FITS image.
    """
    if not isinstance(image, Fits_Simple):
        image = Fits_Simple(image)

    data_masked = image.masked_array
    fig = pyplot.figure(figsize=(8, 6))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_title(image)

    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(data_masked)
    cmap = pyplot.get_cmap()
    cmap.set_bad('r', alpha=0.5)
    ax.imshow(data_masked, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
    pyplot.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax)
    pyplot.show()


def display_many_nickel(path_list):
    """
    Display the data of all images in a list of directories or files.

    Parameters
    ----------
    path_list : list of str
        A list of directory or file paths containing FITS images.
    """
    images = unzip_directories(path_list, output_format='Fits_Simple')
    for image in images:
        display_nickel(image)


def show_image(image: Union[str, Path, CCDData]):
    """
    Display the data of a FITS image using zscale coloring, with masked
    pixels colored red.

    Parameters
    ----------
    image : str, :class:`pathlib.Path`, :class:`~pipnick.utils.fits_class.Fits_Simple`
        The Fits_Simple object or path to the FITS image.
    """
    if not isinstance(image, CCDData):
        _image = CCDData.read(image)

    data_masked = image.masked_array
    fig = pyplot.figure(figsize=(8, 6))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_title(image)

    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(data_masked)
    cmap = pyplot.get_cmap()
    cmap.set_bad('r', alpha=0.5)
    ax.imshow(data_masked, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
    pyplot.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax)
    pyplot.show()
   