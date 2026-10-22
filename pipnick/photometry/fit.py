
###################################################################################
########  Fit Moffat function PSFs to stamps, plot these stamps and plots  ########
###################################################################################

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

from astropy.io import fits
from astropy.modeling.functional_models import Moffat1D
from astropy.visualization import AsinhStretch, ZScaleInterval, ImageNormalize
from astropy.stats import SigmaClip
from astropy.nddata import CCDData
from photutils.aperture import CircularAperture

from pipnick import log
from pipnick.utils.fits_class import Fits_Simple
from pipnick.photometry.moffat import FitEllipticalMoffat2D, FitMoffat2D, make_ellipse
#from pipnick.photometry.starfind import generate_stamps
from pipnick.utils.fits_class import unzip_directories


def fit_psf_stack(input_base, num_images, fittype='ellip', ofile=None):
    """
    Fit one PSF to the stack of all sources found in the directory specified,
    and save this information to relevant files.

    Parameters
    ----------
    input_base : :class:`pathlib.Path`
        Base of path to files with stamp data.
    num_images : int
        Number of images to process.
    fittype : str, optional
        Type of model to fit ('ellip' or 'circ').
    ofile : str, optional
        Output file path.

    Returns
    -------
    :class:`numpy.ndarray`
        Fit of all sources stacked (fit.par = parameters).
    """
    return fit_psf_generic(
        'stack', input_base, num_images, fittype, ofile=ofile
    )

def fit_psf_single(input_base, num_images, fittype='ellip', sigma_clip=True):
    """
    Fit a PSF to each source found in the directory specified, and return the
    source coordinates, fit parameters, and image number.

    Parameters
    ----------
    input_base : :class:`pathlib.Path`
        Base of path to files with stamp data.
    num_images : int
        Number of images to process.
    fittype : str, optional
        Type of model to fit ('ellip' or 'circ').
    sigma_clip : bool, optional
        If True, remove PSF fit sources with unusual FWHM.

    Returns
    -------
    coords : :class:`numpy.ndarray`
        Coordinates of all sources.
    fits : :class:`numpy.ndarray`
        Fits of all sources (fit.par = parameters).
    img_nums : :class:`numpy.ndarray`
        Image number of all sources.
    """
    return fit_psf_generic(
        'single', input_base, num_images, fittype, sigma_clip=sigma_clip, ofile=None
    )

def _moffat_p0(fittype, amplitude, stamp_width, gamma, alpha):
    """
    Construct the initial guess for the Moffat parameters of a stamp.

    Parameters
    ----------
    fittype : str
        Type of model to fit ('ellip' or 'circ').
    amplitude : float
        Initial guess for the amplitude.
    stamp_width : int
        Width of the (square) stamp in pixels; the center is the initial guess
        for the source position.
    gamma : float
        Initial guess for the Moffat core width; both widths for 'ellip'.
    alpha : float
        Initial guess for the Moffat power index.

    Returns
    -------
    :class:`numpy.ndarray`
        The initial parameters, ordered as expected by
        :class:`~pipnick.photometry.moffat.FitEllipticalMoffat2D` ('ellip') or
        :class:`~pipnick.photometry.moffat.FitMoffat2D` ('circ').
    """
    center = float(stamp_width // 2)
    if fittype == 'ellip':
        return np.array([center, center, amplitude, gamma, gamma, 0.0, alpha, 0.0])
    return np.array([center, center, amplitude, gamma, alpha, 0.0])


def fit_psf_generic(
    mode, input_base, num_images, fittype='ellip', sigma_clip=True, ofile=None
):
    """
    Generic function to fit PSFs to images.

    Parameters
    ----------
    mode : str
        Mode of fitting ('stack' or 'single').
    input_base : :class:`pathlib.Path`
        Base of path to files with stamp data.
    num_images : int
        Number of images in directory to process.
    fittype : str, optional
        Type of model to fit ('ellip' or 'circ').
    sigma_clip : bool, optional
        If True, remove PSF fit sources with unusual FWHM.
    ofile : str, optional
        Output file path for 'stack' mode.
    """
    if fittype == 'ellip':
        fitter = FitEllipticalMoffat2D  # Type of fitting function to use
        num_pars = 8
    elif fittype == 'circ':
        fitter = FitMoffat2D  # Type of fitting function to use
        num_pars = 6
    else:
        raise ValueError("fitter must be 'ellip' or 'circ'")

    # Set up directories and files
    # proc_dir = Path('.').resolve() / "proc_files"
    # base = proc_dir / category_str / (category_str) #+ f"_{fittype[:4]}")
    ofits = input_base.with_suffix('.rdx.fits')  # Path to FITS file
    src_ofile = input_base.with_suffix('.src.db')  # Path to source database file

    # Load data from FITS and source files
    hdu = fits.open(ofits)  # Open FITS file
    srcdb = np.genfromtxt(src_ofile, dtype=float)  # Load source database

    # Pull stamp shape and width from FITS file
    stamp_shape = tuple(hdu['STAMPS_01'].data.shape[1:])  # Shape of the stamp images
    stamp_width = stamp_shape[0]  # Width of the stamps
    indx = (srcdb[:,8].astype(int) == 1) & (np.log10(srcdb[:,7]) > 2.0)  # Selection criteria

    # Initialize arrays for 'stack' mode
    if mode == 'stack':
        psf_sum_stack = np.zeros((num_images,) + stamp_shape, dtype=float)
        psf_sum_model = np.zeros((num_images,) + stamp_shape, dtype=float)
        psf_sum_model_par = np.zeros((num_images, num_pars), dtype=float)
    
    # Initialize arrays for 'single' mode
    elif mode == 'single':
        centroid_coords = []  # List to store centroid coordinates
        fit_pars = []  # List to store fit parameters
        fit_objs = []
        source_images = []  # List to store the image number of all sources
        

    # Process each image
    for i in range(num_images):
        on_chip = (srcdb[:,0] == i+1)
        stamp_indx = np.full(on_chip.size, -1, dtype=int)
        stamp_indx[on_chip] = np.arange(np.sum(on_chip))
        ext = f'STAMPS_{i+1:02}'  # Extension name for FITS file
        in_q = on_chip & indx
        
        # Initial parameters for Moffat function
        default_fwhm = 8
        alpha = 3.5
        gamma = default_fwhm / 2 / np.sqrt(2**(1/alpha)-1)
        
        if mode == 'stack':
            # Stack mode: Sum the stamps and divide by flux before fitting
            psf_sum_stack[i,...] = (
                np.sum(hdu[ext].data[stamp_indx[in_q]], axis=0) / np.sum(srcdb[in_q,7])
            )
            # Initial guess for Moffat parameters
            # TODO: The amplitude guess is the maximum of row i of the stack,
            # not of the full stamp, which is likely a bug.
            p0 = _moffat_p0(fittype, np.amax(psf_sum_stack[i][i]), stamp_width, gamma, alpha)
            
            fit = fitter(psf_sum_stack[i])  # Initialize fit object
            fit.fit(p0=p0)  # Perform the fit
            psf_sum_model[i,...] = fit.model()  # Get the model image
            psf_sum_model_par[i,...] = fit.par  # Get the fit parameters
        
        elif mode == 'single':
            # Single mode: Fit each individual stamp
            for step_num, stamp_img in enumerate(hdu[ext].data[stamp_indx[in_q]]):
                # Initial guess for Moffat parameters
                # TODO: The amplitude guess is the maximum of row i of the
                # stamp, not of the full stamp, which is likely a bug.
                p0 = _moffat_p0(fittype, np.amax(stamp_img[i]), stamp_width, gamma, alpha)
                
                fit = fitter(stamp_img)  # Initialize fit object
                try:
                    fit.fit(p0=p0)  # Perform the fit
                except ValueError:
                    continue
                fit_par = fit.par  # Get the fit parameters
                
                # Find centroid coordinates
                condition = (srcdb[:, 0] == i+1) & (srcdb[:, 1] == step_num)
                centroid_x = srcdb[condition][0][2]
                centroid_y = srcdb[condition][0][3]
                
                centroid_coords.append((centroid_x, centroid_y))  # Store coordinates
                fit_pars.append(fit_par)  # Store fit parameters
                fit_objs.append(fit)
                source_images.append(i)
            
    hdu.close()  # Close the FITS file
    # ^ interrupting code during run may leave hdu open--just restart kernel
    
    if mode == 'stack':
        # Save the results to a new FITS file
        fits.HDUList(
            [
                fits.PrimaryHDU(), fits.ImageHDU(data=psf_sum_stack, name='STACK'),
                fits.ImageHDU(data=psf_sum_model, name='MOFFAT'),
                fits.ImageHDU(data=psf_sum_model_par, name='PAR')
            ]
        ).writeto(str(ofile), overwrite=True)
        return fit
    
    elif mode == 'single':
        # Eliminate sources with irregular FWHMs
        fit_pars = np.array(fit_pars)
        fit_objs = np.array(fit_objs)
        source_images = np.array(source_images)
        centroid_coords = np.array(centroid_coords)
        
        if not sigma_clip:
            return centroid_coords, fit_objs, source_images
        else:
            if fittype == 'ellip':
                fwhm1 = FitMoffat2D.to_fwhm(fit_pars[:,3], fit_pars[:,6])
            elif fittype == 'circ':
                fwhm1 = FitMoffat2D.to_fwhm(fit_pars[:,3], fit_pars[:,4])
            
            # Create a SigmaClip object and apply it to get a mask
            sigma_clipper = SigmaClip(sigma=4, maxiters=5)
            masked_fwhm1 = sigma_clipper(fwhm1)
            clipped_fit_pars = fit_pars[~masked_fwhm1.mask]
            clipped_fit_objs = fit_objs[~masked_fwhm1.mask]
            clipped_coords = centroid_coords[~masked_fwhm1.mask]
            clipped_source_images = source_images[~masked_fwhm1.mask]
            
            if fittype == 'ellip':
                fwhm2 = FitMoffat2D.to_fwhm(clipped_fit_pars[:,4], clipped_fit_pars[:,6])
                masked_fwhm2 = sigma_clipper(fwhm2)
                clipped_fit_pars = clipped_fit_pars[~masked_fwhm2.mask]
                clipped_fit_objs = clipped_fit_objs[~masked_fwhm2.mask]
                clipped_coords = clipped_coords[~masked_fwhm2.mask]
                clipped_source_images = clipped_source_images[~masked_fwhm2.mask]
            
            log.info(
                "Number of sources removed in sigma clipping = "
                f"{len(fit_pars) - len(clipped_fit_objs)}"
            )
            log.info(f"Number of sources remaining = {len(clipped_fit_objs)}")
            
            return clipped_coords, clipped_fit_objs, clipped_source_images


#def get_source_pars(path_list, category_str=None, fittype='ellip'):
#    """
#    Extract source coordinates and Moffat fit parameters from image data.
#
#    Parameters
#    ----------
#    path_list : list
#        List of paths (directories or files) to unzip.
#    category_str : str, optional
#        Category string for identifying the path to data.
#    fittype : str, optional
#        Type of model to fit ('ellip' or 'circ').
#
#    Returns
#    -------
#    source_coords : `numpy.ndarray`
#        Array of source coordinates.
#    source_pars : `numpy.ndarray`
#        Array of source parameters (x0, y0, amplitude, gamma1, gamma2, phi, alpha, background).
#    img_nums : `numpy.ndarray`
#        Array of image numbers corresponding to the sources.
#    """
#    # Unzip directories to get image files
#    images = unzip_directories(path_list, output_format='Path')
#
#    # Create output directories
#    proc_dir = Path('.').resolve() / "proc_files"
#    Path.mkdir(proc_dir, exist_ok=True)
#    proc_subdir = proc_dir / fittype
#    Path.mkdir(proc_subdir, exist_ok=True)
#    base_parent = proc_subdir / category_str
#    Path.mkdir(base_parent, exist_ok=True)
#    base = proc_subdir / category_str / category_str
#
#    # Generate stamps (image of sources) for image data
#    generate_stamps(images, output_base=base)
#
#    # Fit PSF models and get source coordinates and parameters
#    source_coords, source_fits, img_nums = fit_psf_single(base, len(images))
#    source_pars = np.array([fit.par for fit in source_fits])
#    return source_coords, source_pars, img_nums


def get_graphable_pars(file_paths, group_name, verbose=False):
    """
    Fit PSF and extract parameters for given files (all stars in these files are
    stacked), storing intermediates in folder proc_files/elliptical/group_name.

    Parameters
    ----------
    file_paths : list
        List of file paths to analyze.
    group_name : str
        Folder for storing intermediates (proc_files/group_name).
    verbose : bool, optional
        If True, print detailed output during processing.

    Returns
    -------
    fwhm : float
        FWHM.
    ecc : float
        FWHM eccentricity.
    phi : float
        Rotation angle phi.
    """
    # Define directory and base path for processed files
    proc_dir = Path('.').resolve() / "proc_files"
    Path.mkdir(proc_dir, exist_ok=True)
    proc_subdir = proc_dir / 'ellip'
    Path.mkdir(proc_subdir, exist_ok=True)
    base_parent = proc_subdir / group_name
    Path.mkdir(base_parent, exist_ok=True)
    base = proc_subdir / group_name / group_name

    # Generate image stamps for the given files
    generate_stamps(file_paths, output_base=base)

    # Fit PSF stack and get the fit results
    psf_file = Path(f'{str(base)}.psf.fits').resolve()  # PSF info stored here
    fit = fit_psf_stack(base, 1, fittype='ellip', ofile=psf_file)

    # Plot PSF and get FWHM and phi values
    plot_file = Path(f'{str(base)}.psf.pdf').resolve()  # Plots stored here
    psf_plot(plot_file, fit, fittype='ellip', plot_fit=True)
    fwhm = get_param_list('fwhm', np.array([fit.par]), (1,))[0][0]
    ecc = get_param_list('ecc', np.array([fit.par]), (1,))[0][0]
    phi = get_param_list('phi', np.array([fit.par]), (1,))[0][0]

    if verbose:
        print(f"Avg FWHM = {fwhm:3f}")
        print(f"FWHM_ecc = {ecc:3f}")
        print(f"Rotation angle phi = {phi:3f}")

    return fwhm, ecc, phi


def get_param_list(param_type, pars, shape, img_nums=None, plate_scale=None):
    """
    Generate the desired single parameter list, color range, and title
    for contour plotting based on Moffat pars.

    Parameters
    ----------
    param_type : str
        Type of parameter ('fwhm', 'fwhm residuals', 'phi', 'ecc').
    pars : :class:`numpy.ndarray`
        Fit parameters (list of par).
    shape : tuple
        Shape to output param_list.
    img_nums : :class:`numpy.ndarray`, optional
        Image number for each source.
    plate_scale : float, optional
        Plate scale of the images in arcsec per *binned* pixel, used to convert
        the FWHM to arcsec.  This is the camera plate scale (see
        :attr:`~pipnick.cameras.cameras.NickelCamera.plate_scale`) times the
        binning.  If None, the FWHM is in pixels.

    Returns
    -------
    param_list : :class:`numpy.ndarray`
        List of parameter values.
    color_range : list
        Range of colors for plotting.
    title : str
        Title for the plot.

    Raises
    ------
    :exc:`ValueError`
        If the input `param_type` is not 'fwhm', 'phi', 'ecc', or 'fwhm residuals'.
    """
    # The FWHM color ranges were optimized for the 2024-06-26 data from the
    # original Nickel camera, which have 0.37 arcsec pixels
    ref_scale = 0.37
    scale, units = (1.0, 'pixels') if plate_scale is None else (plate_scale, 'arcsec')
    if param_type == 'fwhm':
        # Calculate FWHM (average between semi-major and minor axes)
        param_list = (
            FitMoffat2D.to_fwhm(pars[:, 3], pars[:, 6]) + FitMoffat2D.to_fwhm(
                pars[:, 4], pars[:, 6]
            )
        ) / 2 * scale
        color_range = [1.5 / ref_scale * scale, 2.7 / ref_scale * scale]
        title = f"FWHM ({units})"
    elif param_type == 'fwhm residuals':
        # Calculate FWHM residual (relative to minimum FWHM in image)
        fwhm_list = (
            FitMoffat2D.to_fwhm(pars[:, 3], pars[:, 6]) + FitMoffat2D.to_fwhm(
                pars[:, 4], pars[:, 6]
            )
        ) / 2
        mins = {
            img_num: np.min(fwhm_list[img_nums == img_num]) for img_num in list(set(img_nums))
        }
        param_list = np.array(
            [fwhm_list[i] - mins[img_num] for i, img_num in enumerate(img_nums)]
        ) * scale
        color_range = [0.0, 0.36 / ref_scale * scale]
        title = f"FWHM Residuals ({units})"
    elif param_type == 'phi':
        # Convert phi rotation angle relative to x-axis from the original phi
        param_list = np.array(
            [FitEllipticalMoffat2D.get_nice_phi(smooth_par) for smooth_par in pars]
        )
        color_range = [-45., 45.]
        title = "Phi Rotation Angle (deg)"
    elif param_type == 'ecc':
        # Calculate eccentricity
        param_list = []
        for smooth_par in pars:
            fwhm1 = FitMoffat2D.to_fwhm(smooth_par[3], smooth_par[6])
            fwhm2 = FitMoffat2D.to_fwhm(smooth_par[4], smooth_par[6])
            param_list.append(np.sqrt(np.abs(fwhm1 ** 2 - fwhm2 ** 2)) / max(fwhm1, fwhm2))
        param_list = np.array(param_list)
        color_range = [0.29, 0.65]  # Optimized for Nickel 06-26-24 data
        title = "Eccentricity"
    else:
        raise ValueError("Input param_type must be 'fwhm' or 'phi'")

    param_list = param_list.reshape(shape)
    return param_list, color_range, title


def psf_plot(plot_file, fit, fittype='ellip', show=False, plot_fit=True):
    """
    Plot the PSF fitting results and save to a PDF.

    Parameters
    ----------
    plot_file : :class:`pathlib.Path`
        Output PDF file path.
    fit : :class:`~.FitMoffat2D`, :class:`~.FitEllipticalMoffat2D`
        Fitting results.
    fittype : str, optional
        Type of model fit.  Currently only 'ellip' is implemented.
    show : bool, optional
        Show the plot.
    plot_fit : bool, optional
        Include the model fit in the plot.
    """
    if fittype != 'ellip':
        raise ValueError(f"psf_plot() not yet implemented for fittype={fittype}")
    with PdfPages(plot_file) as pdf:
        # Set up the figure
        w, h = plt.figaspect(1.)
        fig = plt.figure(figsize=(1.5*w,1.5*h))
        plt.suptitle(plot_file.stem)  # Set the title of the plot

        stack = fit.c  # Observed stack
        model = fit.model()  # Model stack

        amp = fit.par[2]  # Amplitude of the fit
        if isinstance(fit, FitMoffat2D):
            beta = fit.par[4]  # Moffat beta parameter
            fwhm1 = FitMoffat2D.to_fwhm(fit.par[3], beta)  # Calculate FWHM
            ell_x, ell_y = make_ellipse(fwhm1, fwhm1, 0.)  # Create ellipse for plotting
        else:
            beta = fit.par[6]  # Moffat beta parameter
            phi = fit.get_nice_phi(fit.par)  # Calculate rotation angle
            fwhm1 = FitMoffat2D.to_fwhm(fit.par[3], beta)  # Calculate FWHM1
            fwhm2 = FitMoffat2D.to_fwhm(fit.par[4], beta)  # Calculate FWHM2
            ell_x, ell_y = make_ellipse(fwhm1, fwhm2, fit.par[5])  # Create ellipse for plotting
            if fwhm1 < fwhm2:
                fwhm1, fwhm2 = fwhm2, fwhm1  # Swap FWHMs if necessary
        ell_x += fit.par[0]  # Offset ellipse in x
        ell_y += fit.par[1]  # Offset ellipse in y
            
        # Normalize the images for better visualization
        norm = ImageNormalize(
            np.concatenate((stack, model, stack-model)), interval=ZScaleInterval(contrast=0.10),
            stretch=AsinhStretch()
        )

        # Plot observed stack
        ax = fig.add_axes([0.03, 0.7, 0.2, 0.2])
        ax.imshow(stack, origin='lower', interpolation='nearest', norm=norm)
        ax.contour(stack, [amp/8, amp/4, amp/2, amp/1.1], colors='k', linewidths=0.5)
        ax.set_axis_off()
        ax.text(0.5, 1.01, 'Observed', ha='center', va='bottom', transform=ax.transAxes)

        # Plot model stack
        ax = fig.add_axes([0.24, 0.7, 0.2, 0.2])
        ax.imshow(model, origin='lower', interpolation='nearest', norm=norm)
        ax.contour(model, [amp/8, amp/4, amp/2, amp/1.1], colors='k', linewidths=0.5)
        ax.plot(ell_x, ell_y, color='C3', lw=0.5)
        ax.set_axis_off()
        ax.text(0.5, 1.01, 'Model', ha='center', va='bottom', transform=ax.transAxes)

        # Plot residuals
        ax = fig.add_axes([0.45, 0.7, 0.2, 0.2])
        ax.imshow(stack-model, origin='lower', interpolation='nearest', norm=norm)
        ax.contour(stack-model, [-amp/40, amp/40], colors=['w','k'], linewidths=0.5)
        ax.set_axis_off()
        ax.text(0.5, 1.01, 'Residual', ha='center', va='bottom', transform=ax.transAxes)

        # Plot 1D profiles of source data and model
        r = np.sqrt((fit.x - fit.par[0])**2 + (fit.y - fit.par[1])**2).ravel()
        rlim = np.array([0, 5*fwhm1])  # Radius limits
    
        oned = Moffat1D()  # Initialize 1D Moffat function
        r_mod = np.linspace(*rlim, 100)  # Radial positions
        if plot_fit:
            if isinstance(fit, FitMoffat2D):
                models = [oned.evaluate(r_mod, amp, 0., fit.par[3], beta) + fit.par[5]]
            else:
                models = [
                    oned.evaluate(r_mod, amp, 0., fit.par[3], beta) + fit.par[7],
                    oned.evaluate(r_mod, amp, 0., fit.par[4], beta) + fit.par[7]
                ]
        
        ax = fig.add_axes([0.66, 0.7, 0.3, 0.2])
        ax.minorticks_on()
        ax.set_xlim(rlim)
        ax.tick_params(axis='x', which='both', direction='in')
        ax.tick_params(axis='y', which='both', left=False, right=False)
        ax.yaxis.set_major_formatter(ticker.NullFormatter())
        if plot_fit:
            for model in models:
                ax.plot(r_mod, model, color='C3')
        ax.scatter(r, stack.ravel(), marker='.', lw=0, s=30, alpha=0.5, color='k')

        if isinstance(fit, FitMoffat2D):
            ax.text(
                0.95, 0.9, f'FWHM = {fwhm1:.1f} pix', ha='right', va='center',
                transform=ax.transAxes
            )
            ax.text(
                0.95, 0.78, f'beta = {beta:.2f}', ha='right', va='center', transform=ax.transAxes
            )
        else:
            ax.text(
                0.95, 0.9, f'FWHM_1 = {fwhm1:.1f} pix', ha='right', va='center',
                transform=ax.transAxes
            )
            ax.text(
                0.95, 0.78, f'FWHM_2 = {fwhm2:.1f} pix', ha='right', va='center',
                transform=ax.transAxes
            )
            ax.text(
                0.95, 0.66, f'beta = {beta:.2f}', ha='right', va='center', transform=ax.transAxes
            )
        ax.text(0.5, -0.15, 'R [pix]', ha='center', va='top', transform=ax.transAxes)
            
        pdf.savefig()  # Save the figure to the PDF
        if show:
            plt.show()  # Display the plot
            try:
                fig.clear()
            except:
                pass
        plt.close()


def get_apertures(phot_data, fwhm, scale):
    """Create circular apertures for the sources based on their positions."""
    x = phot_data['x_fit']
    y = phot_data['y_fit']
    return CircularAperture(np.transpose((x, y)), r=2 * fwhm * scale)


def plot_sources(
    image, phot_data, fwhm, flux_name='flux_fit', x_name='x_fit', y_name='y_fit',
    label_name='group_id', scale=1
):
    """
    Plots sources from a photometric data table on a corresponding image, highlighting 
    the grouped/ungrouped sources.

    Parameters
    ----------
    image : str, :class:`astropy.nddata.CCDData`
        Image to plot the sources on. If None, it is loaded from the metadata.
    phot_data : :class:`astropy.table.Table`
        Photometric data containing positions and fluxes of sources.
    given_fwhm : float
        Full-width half-maximum (FWHM) of the sources to set aperture sizes.
    flux_name : str, optional
        Name of the flux column in phot_data (default is 'flux_fit').
    x_name : str, optional
        Name of the x-coordinate column in phot_data (default is 'x_fit').
    y_name : str, optional
        Name of the y-coordinate column in phot_data (default is 'y_fit').
    label_name : str, optional
        Name of the label column in phot_data to use for annotating sources
        (default is 'group_id').
    scale : float, optional
        Scaling factor for the aperture sizes and annotation text (default is 1).
    """
    _image = image if isinstance(image, CCDData) else CCDData.read(image)    
#    if flux_name == 'flux_fit' and 'flux_fit' not in phot_data.colnames:
#        flux_name = 'flux_psf'
    
    # Separate good and bad photometric data based on group size
    is_single = phot_data['group_size'] <= 1
    nsingle = np.sum(is_single)
    in_group = phot_data['group_size'] > 1
    single_apertures = get_apertures(phot_data[is_single], fwhm, scale)
    group_apertures = get_apertures(phot_data[in_group], fwhm, scale)
    
    # Determine image display limits using ZScale
    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(_image)
    
    # Set colormap and mask bad pixels with red
    cmap = plt.get_cmap()
    cmap.set_bad('r', alpha=0.5)
    
    # Plot the image and the good/bad sources
    plt.figure(figsize=(12,10))
#    plt.title(image)
    plt.imshow(
        _image, origin='lower', vmin=vmin, vmax=vmax, cmap=cmap, interpolation='nearest'
    )
    plt.colorbar()
    single_apertures.plot(color='purple', lw=1.5*scale, alpha=1)
    group_apertures.plot(color='r', lw=1.5*scale, alpha=1)
    
    # Annotate singular sources with label_name and flux_name values
    y_offset = 3.5 * fwhm * scale
    _phot_data = phot_data[is_single]
    for i in range(nsingle):
        plt.text(
            _phot_data[x_name][i], _phot_data[y_name][i]+y_offset,
            f'{_phot_data[label_name][i]}: {_phot_data[flux_name][i]:.0f}', color='white',
            fontsize=8*scale, ha='center', va='center'
        )
    
    # Annotate grouped sources with label_name and flux_name values in one large stack
    _phot_data = phot_data[in_group]
    group_ids = np.unique(_phot_data[label_name].data)
    for id in group_ids:
        group = _phot_data[_phot_data[label_name] == id]
        group_x = np.median(group[x_name])
        group_y = np.median(group[y_name]) + y_offset
        for i in range(len(group)):
            plt.text(
                group_x, group_y+i*20*scale, f'{id}: {group[flux_name][i]:.0f}', color='red',
                fontsize=8*scale, ha='center', va='center'
            )
    
    # Show plot
#    plt.gcf().set_dpi(300)
    plt.show()
