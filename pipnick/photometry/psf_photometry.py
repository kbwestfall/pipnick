from pathlib import Path

from IPython import embed

import numpy as np
from scipy.spatial import KDTree

from astropy.table import Table
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.functional_models import Moffat2D
from photutils.detection import IRAFStarFinder
from photutils.psf import IterativePSFPhotometry, make_psf_model
from photutils.background import MMMBackground, MADStdBackgroundRMS, LocalBackground
from photutils.psf import SourceGrouper

from pipnick import logger
from pipnick.utils.fits_class import Fits_Simple
from pipnick.utils.nickel_data import bad_columns
from pipnick.utils.log import log_astropy_table
from pipnick.photometry.moffat_model import MoffatElliptical2D_photutils
from pipnick.photometry import starfind
from pipnick.photometry.fit import fit_psf_single, fit_psf_stack, psf_plot, plot_sources


#np.set_printoptions(edgeitems=100)

def psf_analysis(image, thresh=10.0, mode='all', fittype='circ',
                 plot_final=True, plot_inters=False):
    """
    Perform PSF photometry on a given image, including source detection and PSF fitting.

    This function uses a Moffat function to identify sources in the image,
    generate stamps, fit PSF models, calculate flux, and optionally plot results.

    Parameters
    ----------
    image : `Fits_Simple` or str
        Image for PSF analysis.
    proc_dir : Path
        Directory to store intermediate processing products.
    thresh : float, optional
        Detection threshold in units of background standard deviation. Default is 10.0.
    mode : str, optional
        Mode for PSF fitting in photutils package. Default is 'all'.
    fittype : str, optional
        Type of PSF fit to use ('circ' for circular, 'ellip' for elliptical). Default is 'circ'.
    plot_final : bool, optional
        If True, plot final results. Default is True.
    plot_inters : bool, optional
        If True, plot intermediate results. Default is False.

    Returns
    -------
    `astropy.table.Table`
        Table containing photometric data for detected sources.
    """
    logger.info('Detecting sources')
    source_data = starfind.detect_sources(image, thresh=15., fwhm=7.)

    bkg, stamps = starfind.extract_stamps(image, source_data)

    embed()
    exit()

    # Convert image to Fits_Simple if necessary
    if not isinstance(image, Fits_Simple):
        image = Fits_Simple(image)
    logger.debug(f"analyze_sources() called on image {image.filename}")

    # Prepare a new mask and data, removing bad columns
    new_mask = image.mask.copy().astype(bool)
    new_mask[:5, :] = True  # Mask top 5 rows
    new_mask[-5:, :] = True  # Mask bottom 5 rows
    new_mask[:, :5] = True  # Mask left 5 columns
    new_mask[:, -5:] = True  # Mask right 5 columns
    image.mask = new_mask
    new_data = image.data.copy()
    new_data[:, bad_columns] = 0  # Set bad columns to 0
    image.data = new_data
    img = image.masked_array

    # ----------------------------------------------------------------------
    # Use a Moffat fit to find & fit initial sources
    # ----------------------------------------------------------------------

    # Create output directories
    img_name = image.path.stem.split('_')[0]
    proc_subdir = proc_dir / fittype
    Path.mkdir(proc_subdir, exist_ok=True)  # Create subdirectory for fitting type
    base_parent = proc_subdir / img_name
    Path.mkdir(base_parent, exist_ok=True)  # Create parent directory for image
    base = proc_subdir / img_name / img_name  # Base filename for outputs

    # Generate stamps (image of sources) for image data
    source_data = generate_stamps([image], output_base=base, thresh=thresh)

    # Convert source data into Astropy table
    column_names = ['chip', 'id', 'xcentroid', 'ycentroid', 'bkg', 
                    'kron_radius', 'raw_flux', 'flux', '?']
    sources = Table(source_data, names=column_names)
    table_str = log_astropy_table(sources)
    logger.debug(f"Sources Found (Iter 1): \n{table_str}")

    # Fit PSF models and get source coordinates and parameters
    source_coords, source_fits, _ = fit_psf_single(base, 1, fittype=fittype, sigma_clip=False)
    source_pars = np.array([fit.par for fit in source_fits])

    try:
        # Fit PSF to stack of all sources in image
        psf_file = Path(f'{str(base)}.psf.fits').resolve()  # PSF info stored here
        stack_par = fit_psf_stack(base, 1, fittype=fittype, ofile=psf_file).par
        stack_fwhm = process_par(stack_par, 'Stack', fittype=fittype)
        fit_par = stack_par
        fit_fwhm = stack_fwhm
    except:
        # Handle cases where PSF stack fitting fails -- take avg of first 7 sources' fits
        brightest = np.array(sorted(source_pars, key=lambda coord: coord[2])[:min(7, len(source_pars))])
        clip_avg_par = np.mean(brightest, axis=0)
        clip_avg_fwhm = process_par(clip_avg_par, 'Clipped Avg', fittype=fittype)
        fit_par = clip_avg_par
        fit_fwhm = clip_avg_fwhm

    # Convert initial PSF fit data into Astropy table, calculating flux based on PSF fits
    init_phot_data = Table()
    init_phot_data.add_column(source_coords[:, 0], name='x_fit')
    init_phot_data.add_column(source_coords[:, 1], name='y_fit')
    flux_integrals = [discrete_moffat_integral(par, fittype=fittype, step_size=0.5) for par in source_pars]
    init_phot_data.add_column(flux_integrals, name='flux_fit')
    init_phot_data.add_column(list(range(len(source_pars))), name='group_id')
    init_phot_data.add_column([1] * len(source_pars), name='group_size')
    init_phot_data.meta['image_path'] = image.path

    # Plot intermediate results
    if plot_inters:
        plot_sources(init_phot_data, fit_fwhm)

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
    std = bkgrms(img)

    iraffind = IRAFStarFinder(fwhm=fit_fwhm, threshold=thresh * std, minsep_fwhm=0.1)  # Source finder
    grouper = SourceGrouper(min_separation=2 * fit_fwhm)  # Grouping algorithm
    mmm_bkg = MMMBackground()  # Background estimator
    local_bkg = LocalBackground(*local_bkg_range, mmm_bkg)  # Local background
    fitter = LevMarLSQFitter()  # Optimization algorithm

    # Set bad columns to the background value
    bkgd = mmm_bkg.calc_background(img)
    new_data = img.data.copy()
    new_data[:, bad_columns] = bkgd
    img = np.ma.masked_array(new_data, img.mask)
    image.data = new_data

    # Define the PSF model
    if fittype == 'circ':
        moffat_psf = Moffat2D(gamma=fit_par[3], alpha=fit_par[4])
        moffat_psf = make_psf_model(moffat_psf)
    elif fittype == 'ellip':
        moffat_psf = MoffatElliptical2D_photutils(gamma1=fit_par[3], gamma2=fit_par[4],
                                                  phi=fit_par[5], alpha=fit_par[6])
        moffat_psf = make_psf_model(moffat_psf)

    # Create the photometry object
    phot = IterativePSFPhotometry(finder=iraffind, grouper=grouper,
                                  localbkg_estimator=local_bkg, psf_model=moffat_psf,
                                  fitter=fitter, fit_shape=win,
                                  aperture_radius=aper_size, mode=mode)
    
    # Perform the fitting and collect photometric data
    phot_data = phot(data=img.data, mask=img.mask,
                     init_params=Table(sources['xcentroid', 'ycentroid', 'flux'],
                                       names=('x_0', 'y_0', 'flux_0')))
    phot_data.add_column(image.airmass * np.ones(len(phot_data)), name='airmass')
    phot_data.meta['image_path'] = image.path

    table_str = log_astropy_table(phot_data)
    logger.debug(f"Sources Found (Iter 2): \n{table_str}")

    if plot_inters:
        plot_groups(phot_data, source_coords, source_fits, base)  # Plot group results
    if plot_final:
        plot_sources(phot_data, fit_fwhm)  # Plot final results
    
    return phot_data


def plot_groups(phot_data, source_coords, source_fits, base):
    """
    Show stamp images, side profiles, and PSF fits of sources in groups
    identified in the photometric data.

    Parameters
    ----------
    phot_data : `astropy.table.Table`
        Table containing photometric data of detected sources.
    source_coords : `numpy.ndarray`
        Coordinates of the sources from initial fitting.
    source_fits : list
        List of fitted source models.
    base : Path
        Base path for output plots.
    """
    group_data = phot_data[phot_data['group_size'] > 1]  # Filter groups with multiple sources
    group_ids = list(sorted(set(group_data['group_id'])))  # Unique group IDs

    # Show plots for each group ID
    for id in group_ids:
        logger.warning(f"Group {id} has multiple fitted PSFs: displaying original source")
        group = phot_data[phot_data['group_id'] == id]
        
        # Match coordinates of the group center to original source coordinates
        group_x = np.median(group['x_fit'])
        group_y = np.median(group['y_fit'])
        matching_indices = match_coords((group_x, group_y), source_coords, 2.0)
        if len(matching_indices) == 0:
            matching_indices = match_coords((group_x, group_y), source_coords, 4.0)
            if len(matching_indices) == 0:
                logger.warning("No nearby source found to display")
        if len(matching_indices) > 1:
            logger.info(f"Multiple nearby sources that could match this group; displaying all")
        
        # Display plots for all potential matching sources
        for index in matching_indices:
            matching_fit = source_fits[index]
            plot_file = Path(f'{str(base)}_src{index+1}.psf.pdf').resolve()
            psf_plot(plot_file, matching_fit, show=True, plot_fit=True)


def match_coords(target, search_space, max_dist=2.0):
    """
    Match target coordinates to a search space using KDTree.

    Parameters
    ----------
    target : tuple
        Target coordinates to match.
    search_space : `numpy.ndarray`
        Array of coordinates to search within.
    max_dist : float, optional
        Maximum distance for matching. Default is 2.0.

    Returns
    -------
    list
        List of indices in the search space that are within max_dist of the target.
    """
    search_tree = KDTree(search_space)  # Create KDTree for fast searching
    indices = search_tree.query_ball_point(target, max_dist)  # Find nearby points
    logger.debug(f"Search found indices {indices} within {max_dist} of {target}")
    return indices


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


def discrete_moffat_integral(par, fittype, step_size=1.0):
    """
    Calculate the integral of a Moffat function over a discrete grid.

    Parameters
    ----------
    par : array-like
        Parameters for the Moffat model.
    fittype : str
        Type of Moffat fit ('circ' or 'ellip').
    step_size : float, optional
        Size of the grid step for integration. Default is 1.0.

    Returns
    -------
    float
        Total flux integral over the grid.
    """
    
    # Calculate the start and end points for the grid
    grid_size = 10
    half_size = grid_size // 2
    x_start, x_end = -half_size + step_size / 2, half_size - step_size / 2
    y_start, y_end = half_size - step_size / 2, -half_size + step_size / 2
    x_coords = np.arange(x_start, x_end + step_size, step_size)
    y_coords = np.arange(y_start, y_end - step_size, -step_size)
    grid_x, grid_y = np.meshgrid(x_coords, y_coords)

    # Calculate flux in each pixel based on the fit type
    if fittype == 'circ':
        pixel_fluxes = Moffat2D.evaluate(grid_x, grid_y, par[2], 0, 0, par[3], par[4])
    elif fittype == 'ellip':
        pixel_fluxes = MoffatElliptical2D_photutils.evaluate(grid_x, grid_y, par[2], 0, 0,
                                                             par[3], par[4], par[5], par[6])
    pixel_fluxes *= step_size**2  # Scale by the area of each pixel
    return np.sum(pixel_fluxes)  # Return total flux


def process_par(par, label, fittype):
    """
    Process fit parameters, obtain FWHM, and log relevant information.

    Parameters
    ----------
    par : array-like
        Parameters for the Moffat model.
    label : str
        Label for logging purposes.
    fittype : str
        Type of Moffat fit ('circ' or 'ellip').

    Returns
    -------
    float
        The calculated FWHM based on the fit parameters.
    """
    if fittype == 'circ':
        fwhm = gamma_to_fwhm(par[3], par[4])
        logger.debug(f"{label} Moffat fit parameters: \namplitude = {par[2]}, gamma = {par[3]}, alpha = {par[4]}, bkgd = {par[5]}")
    elif fittype == 'ellip':
        fwhm1 = gamma_to_fwhm(par[3], par[6])
        fwhm2 = gamma_to_fwhm(par[4], par[6])
        fwhm = (fwhm1 + fwhm2) / 2
        logger.debug(f"{label} Moffat fit parameters: \namplitude = {par[2]}, gamma1 = {par[3]}, gamma2 = {par[4]}, phi = {par[5]}, alpha = {par[6]}, bkgd = {par[7]}")
    logger.info(f"{label} FWHM = {fwhm}")
    return fwhm


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
    group_data = group_data.copy()[~np.isin(group_data['group_id'], preserve)]  # Exclude preserved groups

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
    table_str = log_astropy_table(new_data)  # Log consolidated data
    logger.debug(f"Consolidated sources: \n{table_str}")
    return new_data  # Return consolidated data
