"""
Characterize and plot Moffat PSF fits, and plot sources on an image.
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
from matplotlib.figure import Figure

from astropy.modeling.functional_models import Moffat1D
from astropy.visualization import AsinhStretch, ZScaleInterval, ImageNormalize
from astropy.nddata import CCDData
from photutils.aperture import CircularAperture

from pipnick.photometry.moffat import FitEllipticalMoffat2D, FitMoffat2D, make_ellipse


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


def psf_plot(plot_file, fit, plot_fit=True, title=None):
    """
    Plot a PSF fit: the stacked sources, the Moffat model, the residuals, and
    the radial profile.

    The figure is constructed without :mod:`matplotlib.pyplot`, so that it
    can be written in any environment.

    Parameters
    ----------
    plot_file : str, :class:`pathlib.Path`
        Output file; the format is set by its extension (e.g., ``.png`` or
        ``.pdf``).
    fit : object
        The fit, circular or elliptical: a
        :class:`~pipnick.photometry.moffat.FitMoffat2D` or
        :class:`~pipnick.photometry.moffat.FitEllipticalMoffat2D` object.
    plot_fit : bool, optional
        Include the model fit in the radial profile.
    title : str, optional
        The figure title.
    """
    # Set up the figure
    fig = Figure(figsize=(9.6, 3.0))
    if title is not None:
        fig.suptitle(title)

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
    ax = fig.add_axes([0.03, 0.15, 0.2, 0.64])
    ax.imshow(stack, origin='lower', interpolation='nearest', norm=norm)
    ax.contour(stack, [amp/8, amp/4, amp/2, amp/1.1], colors='k', linewidths=0.5)
    ax.set_axis_off()
    ax.text(0.5, 1.01, 'Observed', ha='center', va='bottom', transform=ax.transAxes)

    # Plot model stack
    ax = fig.add_axes([0.24, 0.15, 0.2, 0.64])
    ax.imshow(model, origin='lower', interpolation='nearest', norm=norm)
    ax.contour(model, [amp/8, amp/4, amp/2, amp/1.1], colors='k', linewidths=0.5)
    ax.plot(ell_x, ell_y, color='C3', lw=0.5)
    ax.set_axis_off()
    ax.text(0.5, 1.01, 'Model', ha='center', va='bottom', transform=ax.transAxes)

    # Plot residuals
    ax = fig.add_axes([0.45, 0.15, 0.2, 0.64])
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
    
    ax = fig.add_axes([0.66, 0.15, 0.3, 0.64])
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
    fig.savefig(plot_file)


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
