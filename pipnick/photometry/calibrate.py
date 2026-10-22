import numpy as np
from scipy import linalg

from astropy import stats


def photometric_calibration_coefficients_iterej(
    mag_c, mag_i, airmass, color_defs, order=2, weights=None, gpm=None, rej=3., maxiter=3
):
    r"""
    Iteratively determine the photometric calibration coefficients by
    fitting the data and then rejecting outliers.

    Procedure:

        - Determine the best-fitting coefficients using
          :func:`photometric_calibration_coefficients`.
        - Use :func:`apply_photometric_calibration` to apply the
          coefficents to the instrumental magnitudes.
        - Reject outliers; see ``rej``.
        - Repeat until either no points are rejected or the maximum
          number of iterations is reached.

    Parameters
    ----------
    mag_c : :class:`numpy.ndarray`
        Known catalog stellar magnitudes.  Shape must match ``mag_i``.
    mag_i : :class:`numpy.ndarray`
        Measured instrumental magnitudes, :math:`-2.5\log_{10}(f)` with
        :math:`f` in counts.  Shape is :math:`(N_{\rm obs}, N_{\rm
        band})`, where :math:`N_{\rm obs}` is the number of independent
        observations and :math:`N_{\rm band}` is the number of broadband
        filters used.
    airmass : :class:`numpy.ndarray`
        Airmass of each observation.  Shape must match ``mag_i``.
    color_defs : :class:`numpy.ndarray`
        An array indicating how to construct the color used for
        calibrating each magnitude.  Shape must be :math:`(N_{\rm band},
        2)`.  The array contents provide the column indices in ``mag_c``
        and ``mag_i`` for constructing the color.  E.g.,
        ``color_defs=np.array([[0,1],[2,1], ... ]) means use
        ``mag_c[:,0]-mag_c[:,1]`` to construct the color used to
        calibrate the fluxes in ``mag_i[:,0]``, use
        ``mag_c[:,2]-mag_c[:,1]`` to construct the color used to
        calibrate the fluxes in ``mag_i[:,1]``, etc.
    order : int, optional
        Order of the color terms.  Order=1 adopts only linear terms; see
        :func:`photometric_calibration_coefficients`.
    weights : :class:`numpy.ndarray`, optional
        Weights for each measured magnitude.  If None, weights are
        uniform.  Shape must match ``mag_i``.
    gpm : :class:`numpy.ndarray`, optional
        Boolean flag selecting valid measurements.  Shape must match
        ``mag_i``.
    rej : float, optional
        Number of standard deviations used for rejecting data points
        during each iteration.
    maxiter : int, optional
        Maximum number of fit/rejection iterations.

    Returns
    -------
    coeff : :class:`numpy.ndarray`
        Best-fit calibration coefficients; one set per broadband filter.
    app_mag : :class:`numpy.ndarray`
        Apparent magnitudes constructed by combining the instrumental
        magnitudes with the calibration coefficients.  See
        :func:`apply_photometric_calibration`.  Shape matches ``mag_i``.
    gpm : :class:`numpy.ndarray`
        A boolean array selecting the data used during the calibration.
        Shape matches ``mag_i``.

    """
    # Define the good-point mask
    _gpm = np.ones_like(mag_i.shape, dtype=bool) if gpm is None else gpm
    # ... and its inverse
    mask = np.logical_not(_gpm)
    # ... and the starting number of points reject
    nrej = np.sum(mask)

    # Number of observations and magnitudes
    nobs, nmags = mag_i.shape
    # Instantiate the array used to hold the calibration coefficients
    coeff = np.empty(nmags, dtype=object)

    # Loop until no points are rejected or the maximum number of
    # iterations is reached
    j = 0
    while j < maxiter:
        # Calibrate each band independently
        for i in range(nmags):
            color = np.squeeze(-np.diff(mag_c[:,color_defs[i]], axis=1))
            gpm_i = _gpm[:,i] & _gpm[:,color_defs[i][0]] & _gpm[:,color_defs[i][1]]
            coeff[i] = photometric_calibration_coefficients(
                mag_c[:,i], color, mag_i[:,i], airmass[:,i], order=order, weights=weights,
                gpm=gpm_i
            )
        # Apply the calibration to get the calibrated magnitudes
        app_mag = apply_photometric_calibration(mag_i, airmass, coeff, color_defs, gpm=_gpm)
        # Get the difference between the catalog and calibrated
        # magnitudes, and clip at the requested level
        dmag = stats.sigma_clip(np.ma.MaskedArray(mag_c - app_mag, mask=mask), sigma=rej, axis=0)
        # Update the mask
        mask = np.ma.getmaskarray(dmag)
        # ... determine if any new points were rejected, and break the
        # loop, if not
        if np.sum(mask) - nrej == 0:
            break
        # Update the good-point mask
        _gpm = np.logical_not(mask)
        # ... the number of rejections
        nrej = np.sum(mask)
        # ... and increment the iteration counter
        j += 1

    return coeff, app_mag, _gpm


def photometric_calibration_coefficients(
    mag_c, color_c, mag_i, airmass, order=2, weights=None, gpm=None
):
    r"""
    Determine photometric calibration coefficients.

    Method assumes calibrations of the form

    .. math::

        m - M = k_0 + k_1\ X + k_2\ X C + \sum k_{n+2} C^n

    where :math:`M` are the known apparent magnitudes (``mag_c``),
    :math:`C` are the known colors used for calibration (``color_c``),
    :math:`m` are the measured instrumental magnitudes, :math:`X` is the
    airmass of each observation, and :math:`n\geq1` is the order.

    The approach uses the matrix inversion (linear least-squares) method
    of Harris et al. (1981).

    Parameters
    ----------
    mag_c : :class:`numpy.ndarray`
        Known catalog stellar magnitudes
    color_c : :class:`numpy.ndarray`
        Known catalog stellar colors
    mag_i : :class:`numpy.ndarray`
        Measured instrumental magnitudes, :math:`-2.5\log_{10}(f)` with
        :math:`f` in counts.
    airmass : :class:`numpy.ndarray`
        Airmass of each observation.
    order : int, optional
        Order of the color terms.  Order=1 adopts only linear terms.
        To include quadratic terms, set ``order=2``, etc.
    weights : :class:`numpy.ndarray`, optional
        Weights for each measured magnitude.  If None, weights are
        uniform.
    gpm : :class:`numpy.ndarray`, optional
        Boolean flag selection valid measurements.

    Returns
    -------
    :class:`numpy.ndarray`
        Returns the series of :math:`n+3` coefficients, where :math:`n`
        is the fit ``order``.
    """
    # Setup the boolean mask use to select good measurements
    _gpm = np.ones(mag_i.size, dtype=bool) if gpm is None else gpm
    # Set the number of good measurements
    ngood = np.sum(_gpm)
    # Setup the weights
    _wgt = np.ones(ngood, dtype=float) if weights is None else weights[_gpm]
    # Eqn. 3.2 from Harris et al. (1981) for x_ij
    x = calibration_design_matrix(airmass[_gpm], color_c[_gpm], order)
    # Eqn. 3.10 from Harris et al. (1981) for W_ik
    W = np.dot(_wgt[None,:] * x, x.T)
    # Eqn. 3.9 from Harris et al. (1981) for r_k
    r = np.dot(_wgt[None,:] * x, mag_i[_gpm] - mag_c[_gpm])
    # Eqn. 3.8 from Harris et al. (1981) for a
    return linalg.solve(W, r)


def calibration_design_matrix(airmass, color, order=2):
    """
    Construct the design matrix used for photometric calibration.

    This follows eqn. 3.2 from Harris et al. (1981) for x_ij.

    Parameters
    ----------

    airmass : :class:`numpy.ndarray`
        1D array with the airmass of each observation
    color : :class:`numpy.ndarray`
        1D array with the stellar colors
    order : int, optional
        Order of the color terms.  Order=1 adopts only linear terms.
        To include quadratic terms, set ``order=2``, etc.
   
    Returns
    -------
    x : :class:`numpy.ndarray`
        The result for eqn. 3.2 from Harris et al. (1981) for x_ij.
    """
    x = [np.ones_like(airmass), airmass, airmass*color, color] 
    if order > 1:
        x += [color**i for i in range(2,order+1)]
    return np.stack(x)


def apparent_mag(mag_i, color, airmass, coeffs):
    r"""
    Calculate the apparent magnitudes given the relevant measurements
    and the calibration coefficients.

    This is equivalent to the first equation in the equations sets 2.9
    and 3.1 from Harris et al. (1981).

    Parameters
    ----------
    mag_i : :class:`numpy.ndarray`
        1D vector with the measured instrumental magnitudes,
        :math:`-2.5\log_{10}(f)` with :math:`f` in counts.
    color : :class:`numpy.ndarray`
        Known catalog stellar colors, *not* the instrumental colors.
    airmass : :class:`numpy.ndarray`
        Airmass of each observation.
    coeffs : :class:`numpy.ndarray`
        1D array with the calibration coefficients.  The order of the
        fit is 3 less than the length of this array; i.e., this should
        contain 5 coefficients for an order-2 fit.
    
    Returns
    -------
    :class:`numpy.ndarray`
        1D vector with the calibrated apparent magnitudes.
    """
    x = calibration_design_matrix(airmass, color, len(coeffs)-3)
    return mag_i - np.dot(coeffs, x)


def apply_photometric_calibration(
    mag, airmass, coeff, color_defs, atol=1e-2, maxiter=10, gpm=None
):
    r"""
    Apply photometric calibration coefficients to a set of observations.

    The photometric calibration coefficients are defined by a *known*
    set of stellar colors.  Therefore, this function runs iteratively.
    The first guess at the colors use the instrumental magnitudes.  For
    each iteration, the colors are updated based on the apparent
    magnitudes measured in the previous iteration, and the apparent
    magnitudes are recalculated.  These iterations proceed until the rms
    difference between the original and updated colors is less than
    ``atol`` or the maximum number of iterations (``maxiter``) as been
    reached.

    Parameters
    ----------
    mag : :class:`numpy.ndarray`
        The instrumental magnitudes to calibrate, where
        :math:`-2.5\log_{10}(f)` with :math:`f` in counts.  Shape is
        :math:`(N_{\rm obs}, N_{\rm band})`, where :math:`N_{\rm obs}`
        is the number of independent observations and :math:`N_{\rm
        band}` is the number of broadband filters used.
    airmass : :class:`numpy.ndarray`
        Airmass of each observation.  Shape must match ``mag``.
    coeff : :class:`numpy.ndarray`
        Best-fit calibration coefficients; one set per broadband filter.
    color_defs : :class:`numpy.ndarray`
        An array indicating how to construct the color used for
        calibrating each magnitude.  Shape must be :math:`(N_{\rm band},
        2)`.  The array contents provide the column indices in ``mag_c``
        and ``mag_i`` for constructing the color.  E.g.,
        ``color_defs=np.array([[0,1],[2,1], ... ]) means use
        ``mag_c[:,0]-mag_c[:,1]`` to construct the color used to
        calibrate the fluxes in ``mag_i[:,0]``, use
        ``mag_c[:,2]-mag_c[:,1]`` to construct the color used to
        calibrate the fluxes in ``mag_i[:,1]``, etc.
    atol : float, optional
        The absolute tolerance for the RMS difference between the
        original and updated colors; see function description.
    maxiter : int, optional
        Maximum number of update iterations.
    gpm : :class:`numpy.ndarray`, optional
        Boolean flag selecting valid measurements.  Shape must match
        ``mag``.

    Returns
    -------
    app_mag : :class:`numpy.ndarray`
        Calibrated apparent magnitudes.  Shape matches ``mag``.

    Raises
    ------
    :exc:`ValueError`
        Raised if input arrays have inappropriate/mismatching shapes.

    """
    # Check the input
    if mag.shape != airmass.shape:
        raise ValueError('Shape of mag and airmass arrays must match.')
    nobs, nmag = mag.shape
    if len(coeff) != nmag:
        raise ValueError('Must provide one set of coefficients for each magniude.')
    if color_defs.shape != (nmag, 2):
        raise ValueError('Must define a color used by each magnitude.')

    # Calculate the colors
    color = np.column_stack([-np.diff(mag[:,c], axis=1) for c in color_defs])

    # Ensure all magnitudes are valid for each measurement
    # TODO: Actually only need to make sure each color (magnitude pair) is valid...
    if gpm is not None:
        _gpm = np.all(gpm, axis=1)

    # Iteratively calculate the apparent mags until the colors are as close as specified.
    rms = atol + 1.
    i = 0
    while rms > atol and i < maxiter:
        # Calculate all the apparent mags
        app_mag = np.column_stack(
            [apparent_mag(m, c, a, k) for m, c, a, k in zip(mag.T, color.T, airmass.T, coeff)]
        )
        # Get the updated color
        updated_color = np.column_stack([-np.diff(app_mag[:,c], axis=1) for c in color_defs])
        # Get the RMS difference with the old color
        if gpm is None:
            rms = np.sqrt(np.mean((color-updated_color)**2))
        else:
            rms = np.sqrt(np.mean((color[_gpm]-updated_color[_gpm])**2))
        # Save the updated value
        color = updated_color
        # Increment the iteration counter
        i += 1
    # Return the apparent magnitudes
    return app_mag


