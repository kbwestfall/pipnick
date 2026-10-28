"""
Read the images written by pipnick.

The reduced images (see :func:`~pipnick.pipelines.reduction.reduce_all`) and
the mosaics (see :func:`~pipnick.pipelines.mosaic.mosaic_all`) have the image
in the primary extension, a ``MASK`` extension (non-zero for pixels to
ignore), and an ``UNCERT`` extension with the error of each pixel (a standard
deviation).
"""

import numpy as np

from astropy.io import fits

from pipnick.pkg.exceptions import PipnickDataModelError


def read_reduced(path):
    """
    Read an image, its mask, its error, and its header.

    The data are returned as they are in the file; e.g., the masked pixels
    are not changed.

    Parameters
    ----------
    path : str, :class:`pathlib.Path`
        The FITS file.

    Returns
    -------
    data : :class:`numpy.ndarray`
        The image (``float``).
    mask : :class:`numpy.ndarray`
        Boolean array selecting the masked pixels (the ``MASK`` extension), or
        no pixels if the file has no mask.
    error : :class:`numpy.ndarray`
        The error of each pixel (the ``UNCERT`` extension), or None if the
        file has no error.
    header : :class:`astropy.io.fits.Header`
        The primary header.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the error is not a standard deviation.
    """
    with fits.open(path, memmap=False) as hdu:
        header = hdu[0].header.copy()
        data = np.asarray(hdu[0].data, dtype=float)
        if 'MASK' in hdu:
            mask = np.asarray(hdu['MASK'].data) != 0
        else:
            mask = np.zeros(data.shape, dtype=bool)
        error = None
        if 'UNCERT' in hdu:
            utype = hdu['UNCERT'].header.get('UTYPE', 'StdDevUncertainty')
            if utype != 'StdDevUncertainty':
                raise PipnickDataModelError(
                    f'The UNCERT extension of {path} is a {utype}, not a standard deviation.'
                )
            error = np.asarray(hdu['UNCERT'].data, dtype=float)
    return data, mask, error, header


def invalid_error(error):
    """
    Return the pixels without a valid error.

    Parameters
    ----------
    error : :class:`numpy.ndarray`
        The error of each pixel.

    Returns
    -------
    :class:`numpy.ndarray`
        Boolean array selecting the pixels whose error is not finite or not
        positive; these pixels cannot be weighted, so they are ignored.
    """
    return np.logical_not(np.isfinite(error)) | (error <= 0)
