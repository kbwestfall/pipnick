from pathlib import Path

from IPython import embed

import numpy as np
from scipy.spatial import KDTree

from astropy.table import Table
from astropy.nddata import CCDData
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.functional_models import Moffat2D
from photutils.detection import IRAFStarFinder
from photutils.psf import IterativePSFPhotometry, make_psf_model, MoffatPSF
from photutils.background import MMMBackground, MADStdBackgroundRMS, LocalBackground
from photutils.psf import SourceGrouper

from pipnick import logger
from pipnick.utils.fits_class import Fits_Simple
from pipnick.utils.nickel_data import bad_columns
from pipnick.utils.log import log_astropy_table
from pipnick.photometry.moffat import FitMoffat2D, FitEllipticalMoffat2D, EllipticalMoffat2DModel
from pipnick.photometry import starfind
from pipnick.photometry.fit import fit_psf_single, fit_psf_stack, psf_plot, plot_sources


#np.set_printoptions(edgeitems=100)

def build_moffat_psf_model(image, thresh=15.0, fwhm=7., fittype='circ', pdf_file=None):
    """
    Use well-detected and well-separated sources to model the image PSF as a
    2D Moffat function.

    Parameters
    ----------
    image : str, CCDData
        Image for PSF analysis.
    thresh : float, optional
        Detection threshold in units of background standard deviation.
    fwhm : float, optional
        Initial estimate of the FWHM
    mode : str, optional
        Mode for PSF fitting in photutils package.
    fittype : str, optional
        Type of PSF fit to use ('circ' for circular, 'ellip' for elliptical).

    Returns
    -------
    stack : numpy.ma.MaskedArray
        The stacked source image
    par : numpy.ndarray
        The best-fit model parameters
    model : numpy.ndarray
        The best-fit model sampled at each pixel of the image stamps
    """
    # Read the image if necessary
    _image = image if isinstance(image, CCDData) else CCDData.read(image)
    source_data = starfind.detect_sources(_image, thresh=thresh, fwhm=fwhm)
    logger.info(f'Found {len(source_data)} well-separated sources with well-detected fluxes.')

    logger.info('Extracting image stamps.')
    bkg, stamps = starfind.extract_stamps(_image, source_data)
    npix_stamp = np.prod(stamps.shape[1:])
    mask_fraction = np.sum(np.ma.getmaskarray(stamps), axis=(1,2)) / npix_stamp
    indx = source_data['good'].data & (mask_fraction < 0.01)
    logger.info(f'Found {np.sum(indx)} good sources to stack for the PSF measurement.')

    if pdf_file is not None:
        starfind.plot_stamps(stamps[indx], pdf_file, vmin=None, vmax=None)

    # TODO: mean or median?
    stack = np.ma.mean(stamps[indx], axis=0)

    if fittype == 'circ':
        fitter = FitMoffat2D(stack) 
    elif fittype == 'ellip':
        fitter = FitEllipticalMoffat2D(stack)
    else:
        logger.error(f'Unknown fit type {fittype}')

    fitter.fit()
    return stack, fitter.par, fitter.model()


def psf_analysis(image, thresh=10.0, guess_fwhm=7., mode='all', fittype='circ', show=False):
    """
    Perform PSF photometry on a given image, including source detection and PSF fitting.

    This function uses a Moffat function to identify sources in the image,
    generate stamps, fit PSF models, calculate flux, and optionally plot results.

    Parameters
    ----------
    image : str, CCDData
        Image for PSF analysis.
    thresh : float, optional
        Detection threshold in units of background standard deviation.
    guess_fwhm : float, optional
        Initial guess for the FWHM of the PSF in pixels
    mode : str, optional
        Mode for PSF fitting in photutils package. Default is 'all'.
    fittype : str, optional
        Type of PSF fit to use ('circ' for circular, 'ellip' for elliptical). Default is 'circ'.
    show : bool, optional
        Show the image with the sources marked

    Returns
    -------
    numpy.ndarray
        Best-fit parameters for the PSF used for the photometry
    `astropy.table.Table`
        Table containing the results of the PSF photometry
    """
    _image = image if isinstance(image, CCDData) else CCDData.read(image)
    psf_stack, moffat_par, moffat_model = build_moffat_psf_model(
        _image, thresh=thresh, fwhm=guess_fwhm, fittype=fittype
    )
    fit_fwhm = process_par(moffat_par, fittype=fittype)
    logger.info(f'Measured Moffat PSF FWHM = {fit_fwhm:.1f} pixels')

    # ----------------------------------------------------------------------
    # Attempt to improve the source detection by improving the FWHM estimate
    # ----------------------------------------------------------------------
    aper_size = fit_fwhm * 1.8  # Set aperture size based on FWHM
    local_bkg_range = (3 * fit_fwhm, 6 * fit_fwhm)  # Update background range
    win = int(np.ceil(2 * fit_fwhm))  # Calculate window size for fitting
    if win % 2 == 0:
        win += 1  # Ensure window size is odd

    # Background estimation and source finding
    bkgrms = MADStdBackgroundRMS()
    std = bkgrms(_image)

    iraffind = IRAFStarFinder(fwhm=fit_fwhm, threshold=thresh * std, minsep_fwhm=0.1)  # Source finder
    # TODO: Make this a parameter
    grouper = SourceGrouper(min_separation=2 * fit_fwhm)  # Grouping algorithm
    mmm_bkg = MMMBackground()  # Background estimator
    local_bkg = LocalBackground(*local_bkg_range, bkg_estimator=mmm_bkg)  # Local background
    fitter = LevMarLSQFitter()  # Optimization algorithm

    # Define the PSF model
    if fittype == 'circ':
        moffat_psf = MoffatPSF(alpha=moffat_par[3], beta=moffat_par[4])
#        moffat_psf = Moffat2D(gamma=moffat_par[3], alpha=moffat_par[4])
#        moffat_psf = make_psf_model(moffat_psf)
    elif fittype == 'ellip':
        # TODO: Make this into a "PSF" model similar to "MoffatPSF"
        moffat_psf = EllipticalMoffat2DModel(
            gamma1=moffat_par[3], gamma2=moffat_par[4], phi=moffat_par[5], alpha=moffat_par[6]
        )
        moffat_psf = make_psf_model(moffat_psf)

    # Create the photometry object
    phot = IterativePSFPhotometry(finder=iraffind, grouper=grouper,
                                  localbkg_estimator=local_bkg, psf_model=moffat_psf,
                                  fitter=fitter, fit_shape=win,
                                  aperture_radius=aper_size, mode=mode)
    
    # Perform the fitting and collect photometric data
    phot_data = phot(data=_image.data, mask=_image.mask)
    phot_data['airmass'] = np.full(len(phot_data), _image.header['AIRMASS'], dtype=float)

    if show:
        plot_sources(_image, phot_data, fit_fwhm)

    return moffat_par, phot_data


# def plot_groups(phot_data, source_coords, source_fits, base):
#     """
#     Show stamp images, side profiles, and PSF fits of sources in groups
#     identified in the photometric data.
# 
#     Parameters
#     ----------
#     phot_data : `astropy.table.Table`
#         Table containing photometric data of detected sources.
#     source_coords : `numpy.ndarray`
#         Coordinates of the sources from initial fitting.
#     source_fits : list
#         List of fitted source models.
#     base : Path
#         Base path for output plots.
#     """
#     group_data = phot_data[phot_data['group_size'] > 1]  # Filter groups with multiple sources
#     group_ids = list(sorted(set(group_data['group_id'])))  # Unique group IDs
# 
#     # Show plots for each group ID
#     for id in group_ids:
#         logger.warning(f"Group {id} has multiple fitted PSFs: displaying original source")
#         group = phot_data[phot_data['group_id'] == id]
#         
#         # Match coordinates of the group center to original source coordinates
#         group_x = np.median(group['x_fit'])
#         group_y = np.median(group['y_fit'])
#         matching_indices = match_coords((group_x, group_y), source_coords, 2.0)
#         if len(matching_indices) == 0:
#             matching_indices = match_coords((group_x, group_y), source_coords, 4.0)
#             if len(matching_indices) == 0:
#                 logger.warning("No nearby source found to display")
#         if len(matching_indices) > 1:
#             logger.info(f"Multiple nearby sources that could match this group; displaying all")
#         
#         # Display plots for all potential matching sources
#         for index in matching_indices:
#             matching_fit = source_fits[index]
#             plot_file = Path(f'{str(base)}_src{index+1}.psf.pdf').resolve()
#             psf_plot(plot_file, matching_fit, show=True, plot_fit=True)
# 
# 
# def match_coords(target, search_space, max_dist=2.0):
#     """
#     Match target coordinates to a search space using KDTree.
# 
#     Parameters
#     ----------
#     target : tuple
#         Target coordinates to match.
#     search_space : `numpy.ndarray`
#         Array of coordinates to search within.
#     max_dist : float, optional
#         Maximum distance for matching. Default is 2.0.
# 
#     Returns
#     -------
#     list
#         List of indices in the search space that are within max_dist of the target.
#     """
#     search_tree = KDTree(search_space)  # Create KDTree for fast searching
#     indices = search_tree.query_ball_point(target, max_dist)  # Find nearby points
#     logger.debug(f"Search found indices {indices} within {max_dist} of {target}")
#     return indices


def gamma_to_fwhm(gamma, alpha):
    """
    Convert gamma to full-width half-maximum (FWHM).

    Parameters
    ----------
    gamma : float
        Moffat parameter gamma.
    alpha : float
        Moffat parameter alpha.

    Returns
    -------
    float
        FWHM value calculated from gamma and alpha.
    """
    return 2 * gamma * np.sqrt(2**(1/alpha) - 1)


#def discrete_moffat_integral(par, fittype, step_size=1.0):
#    """
#    Calculate the integral of a Moffat function over a discrete grid.
#
#    Parameters
#    ----------
#    par : array-like
#        Parameters for the Moffat model.
#    fittype : str
#        Type of Moffat fit ('circ' or 'ellip').
#    step_size : float, optional
#        Size of the grid step for integration. Default is 1.0.
#
#    Returns
#    -------
#    float
#        Total flux integral over the grid.
#    """
#    
#    # Calculate the start and end points for the grid
#    grid_size = 10
#    half_size = grid_size // 2
#    x_start, x_end = -half_size + step_size / 2, half_size - step_size / 2
#    y_start, y_end = half_size - step_size / 2, -half_size + step_size / 2
#    x_coords = np.arange(x_start, x_end + step_size, step_size)
#    y_coords = np.arange(y_start, y_end - step_size, -step_size)
#    grid_x, grid_y = np.meshgrid(x_coords, y_coords)
#
#    # Calculate flux in each pixel based on the fit type
#    if fittype == 'circ':
#        pixel_fluxes = Moffat2D.evaluate(grid_x, grid_y, par[2], 0, 0, par[3], par[4])
#    elif fittype == 'ellip':
#        pixel_fluxes = EllipticalMoffat2DModel.evaluate(grid_x, grid_y, par[2], 0, 0,
#                                                             par[3], par[4], par[5], par[6])
#    pixel_fluxes *= step_size**2  # Scale by the area of each pixel
#    return np.sum(pixel_fluxes)  # Return total flux


def process_par(par, fittype):
    """
    Process fit parameters, obtain FWHM, and log relevant information.

    Parameters
    ----------
    par : array-like
        Parameters for the Moffat model.
    fittype : str
        Type of Moffat fit ('circ' or 'ellip').

    Returns
    -------
    float
        The calculated FWHM based on the fit parameters.
    """
    if fittype == 'circ':
        return gamma_to_fwhm(par[3], par[4])
    if fittype == 'ellip':
        return (gamma_to_fwhm(par[3], par[6]) + gamma_to_fwhm(par[4], par[6])) / 2
    logger.error(f'Unknown fit type {fittype}')


def consolidate_groups(phot_data, preserve=[]):
    """
    Consolidate groups of sources in photometric source into one source.

    This function combines data for sources that are identified as belonging 
    to the same group by summing fluxes and weighted-averaging other parameters.

    Parameters
    ----------
    phot_data : `astropy.table.Table`
        Table containing photometric data of detected sources.
    preserve : list, optional
        List of group IDs to preserve during consolidation. Default is an empty list.

    Returns
    -------
    `astropy.table.Table`
        A new table containing consolidated photometric data.
    """
    new_data = phot_data.copy()
    group_data = new_data[new_data['group_size'] > 1]  # Filter groups with multiple sources
    # Exclude preserved groups
    group_data = group_data.copy()[np.logical_not(np.isin(group_data['group_id'], preserve))]

    group_ids = list(sorted(set(group_data['group_id'])))  # Unique group IDs
    logger.info(f"Consolidating groups {group_ids}")
    logger.info(f"Preserving groups {preserve}")

    new_data = phot_data.copy()[~np.isin(phot_data['group_id'], group_ids)]  # Start with non-group data
    
    # Consolidate data for each group
    for id in group_ids:
        group = group_data[group_data['group_id'] == id]
        new_row = {
            'id': id,
            'group_id': id,
            'group_size': group['group_size'][0],
            'iter_detected': group['iter_detected'][-1],
            'local_bkg': np.average(group['local_bkg'], weights=abs(group['flux_fit'])),
            'x_init': np.average(group['x_init'], weights=abs(group['flux_fit'])),
            'y_init': np.average(group['y_init'], weights=abs(group['flux_fit'])),
            'flux_init': np.average(group['flux_init'], weights=abs(group['flux_fit'])),
            'x_fit': np.average(group['x_fit'], weights=abs(group['flux_fit'])),
            'y_fit': np.average(group['y_fit'], weights=abs(group['flux_fit'])),
            'flux_fit': np.sum(group['flux_fit']),
            'x_err': np.sum(group['x_err']),
            'y_err': np.sum(group['y_err']),
            'flux_err': np.sum(group['flux_err']),
            'airmass': group['airmass'][0],
            'npixfit': np.sum(group['npixfit']),
            'qfit': np.average(group['qfit'], weights=abs(group['flux_fit'])),
            'cfit': np.average(group['cfit'], weights=abs(group['flux_fit'])),
            'flags': np.bitwise_or.reduce(group['flags'])  # Combine flags using bitwise OR
        }
        new_data.add_row(new_row)  # Add new row for consolidated group
    
    new_data.sort('group_id')  # Sort consolidated data by group ID
#    table_str = log_astropy_table(new_data)  # Log consolidated data
#    logger.debug(f"Consolidated sources: \n{table_str}")
    return new_data  # Return consolidated data
