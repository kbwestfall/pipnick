import numpy as np

from astropy.nddata import CCDData
from astropy.table import MaskedColumn
from photutils.aperture import CircularAperture
from photutils.aperture import ApertureStats

from pipnick import log
from pipnick.utils.fits_class import Fits_Simple
from pipnick.utils.log import log_astropy_table


def aperture_analysis(image, x, y, aper_size=8.0):
    """
    Perform aperture photometry on an image, based on a PSF photometry
    source catalog.

    This function computes the flux within circular apertures centered 
    at source coordinates, and calculates the ratio of PSF flux to aperture flux.

    Parameters
    ----------
    image : str, :class:`astropy.nddata.CCDData`
        Image from which to perform aperture photometry.
    x, y : :class:`numpy.ndarray`
        Pixel coordinates at which to measure the aperture photometry
    aper_size : float, optional
        Radius of the circular aperture in pixels.

    Returns
    -------
    :class:`astropy.table.Table`
        Table containing updated photometric data, including aperture flux and flux ratios.
    """
    # Ensure image is a Fits_Simple object
    if not isinstance(image, Fits_Simple):
        image = Fits_Simple(image)

    # Sort photometric data by group_id
    phot_data.sort('group_id')

    # Create list of source coordinates for aperture analysis
    coordinates = list(zip(phot_data['x_fit'], phot_data['y_fit']))
    apertures = CircularAperture(coordinates, r=aper_size)
    
    # Calculate aperture statistics including sum within apertures
    aperstats = ApertureStats(image.data, apertures, local_bkg=phot_data['local_bkg'])

    # Copy photometric data for modifications
    all_data = phot_data.copy()
    all_data['flux_fit'].name = 'flux_psf'
    col_index = all_data.colnames.index('flux_psf') + 1
    
    # Add the flux within apertures to the photometric data, masking sources
    # with masked pixels, which aperture photometry cannot handle
    has_masked_pixels = all_data['flags'] % 2 == 1
    all_data.add_column(
        MaskedColumn(aperstats.sum, mask=has_masked_pixels), name='flux_aper', index=col_index
    )
    all_data['flux_aper'].info.format = '.3f'
    log.info(
        "Aperture photometry cannot handle masked pixels--sources with masked pixels "
        "have a masked flux_aper"
    )

    # Calculate the ratio of PSF flux to aperture flux and add it to the table;
    # the ratio is masked where flux_aper is masked
    all_data.add_column(
        all_data['flux_psf'] / all_data['flux_aper'], name='ratio_flux', index=col_index + 1
    )
    all_data['ratio_flux'].info.format = '.3f'
    
    # Log the results of the photometry
    log.debug(f"PSF & Aperture Photometry Results: \n{log_astropy_table(all_data)}")

    return format_table(all_data)


def format_table(phot_data):
    """
    Format the photometric data table to include only relevant columns and set formatting.

    Parameters
    ----------
    phot_data : :class:`astropy.table.Table`
        Table containing photometric data to be formatted.

    Returns
    -------
    :class:`astropy.table.Table`
        Formatted table containing a concise selection of columns.
    """
    # Specify the desired column order and names
    colnames = [
        'group_id', 'group_size', 'flags', 'x_fit', 'y_fit', 'flux_psf', 'flux_aper', 'ratio_flux',
        'local_bkg', 'x_err', 'y_err', 'flux_err', 'airmass', 'id', 'iter_detected', 'npixfit',
        'qfit', 'cfit'
    ]
    
    # Create a concise version of the data table
    concise_data = phot_data[colnames]
    
    # Set formatting for float columns to three decimal places
    for col in colnames:
        if isinstance(concise_data[col][0], np.float64):
            concise_data[col].info.format = '.3f'
    
    return concise_data
