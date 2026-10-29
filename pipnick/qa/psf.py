"""
QA figure of the Moffat model of the PSF of an image.
"""

import numpy as np
from matplotlib import ticker
from matplotlib.figure import Figure

from astropy.modeling.functional_models import Moffat1D
from astropy.visualization import AsinhStretch, ZScaleInterval, ImageNormalize

from pipnick.photometry.moffat import FitMoffat2D, gamma_to_fwhm, make_ellipse
from pipnick.qa.display import save_figure


def plot_psf(ofile, fit, plot_fit=True, title=None):
    """
    Plot a PSF fit: the stacked sources, the Moffat model, the residuals, and
    the radial profile.

    Parameters
    ----------
    ofile : str, :class:`pathlib.Path`
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
        fwhm1 = gamma_to_fwhm(fit.par[3], beta)  # Calculate FWHM
        ell_x, ell_y = make_ellipse(fwhm1, fwhm1, 0.)  # Create ellipse for plotting
    else:
        beta = fit.par[6]  # Moffat beta parameter
        fwhm1 = gamma_to_fwhm(fit.par[3], beta)  # Calculate FWHM1
        fwhm2 = gamma_to_fwhm(fit.par[4], beta)  # Calculate FWHM2
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
    save_figure(fig, ofile)
