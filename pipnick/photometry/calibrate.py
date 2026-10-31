r"""
Photometric calibration.

The calibration of each band relates the instrumental magnitudes, :math:`m`,
to the known (catalog) magnitudes, :math:`M`, of standard stars:

.. math::

    m - M = \sum_t k_t\ x_t,

where the terms :math:`x_t` are selected by name (see
:func:`design_matrix`) from

    - ``zp``: the zero point, :math:`x = 1`;
    - ``airmass``: the extinction, :math:`x = X`;
    - ``color``: the color term, :math:`x = C`;
    - ``color2``: the quadratic color term, :math:`x = C^2`;
    - ``airmass_color``: the color dependence of the extinction,
      :math:`x = X C`;

for the airmass :math:`X` and the (catalog) color :math:`C`.  The
instrumental magnitudes are :math:`m = -2.5\log(f / t)`, for the flux
:math:`f` in electrons and the exposure time :math:`t`, so the zero point
applies to count rates in electrons per second.

The coefficients are found by weighted linear least squares, following the
matrix-inversion method of Harris et al. (1981, PASP, 93, 507).  Any term can
be fixed at a given value, in which case it is subtracted before the fit.
The errors of the coefficients (and those propagated from them) are
statistical only: fixed terms have no error, and the systematic error they
introduce is not included.
"""

from dataclasses import dataclass

import numpy as np

from astropy import stats
from pipnick.pkg.exceptions import PipnickParameterError
from pipnick.utils.lsq import weighted_lstsq

# The terms of the calibration model, in the order used by the design matrix
CALIBRATION_TERMS = ('zp', 'airmass', 'color', 'color2', 'airmass_color')


def check_terms(terms):
    """
    Check and order the terms of a calibration model.

    Parameters
    ----------
    terms : list
        The names of the terms; see :data:`CALIBRATION_TERMS`.  The zero
        point (``zp``) is required.

    Returns
    -------
    list
        The terms, in the order of :data:`CALIBRATION_TERMS`.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if a term is unknown or repeated, or the zero point is missing.
    """
    unknown = [t for t in terms if t not in CALIBRATION_TERMS]
    if len(unknown) > 0:
        raise PipnickParameterError(
            f'Unknown calibration terms: {", ".join(unknown)}.  Options are: '
            f'{", ".join(CALIBRATION_TERMS)}.'
        )
    if len(set(terms)) != len(terms):
        raise PipnickParameterError(f'Calibration terms are repeated: {", ".join(terms)}.')
    if 'zp' not in terms:
        raise PipnickParameterError('The calibration must include the zero point (zp).')
    return [t for t in CALIBRATION_TERMS if t in terms]


def design_matrix(terms, airmass, color):
    """
    Construct the design matrix of a calibration model.

    This follows eqn. 3.2 of Harris et al. (1981) for :math:`x_{ij}`.

    Parameters
    ----------
    terms : list
        The names of the terms; see :func:`check_terms`.
    airmass : :class:`numpy.ndarray`
        The airmass of each measurement.
    color : :class:`numpy.ndarray`
        The color of each measurement.

    Returns
    -------
    :class:`numpy.ndarray`
        The design matrix, with shape (number of terms, number of
        measurements).
    """
    _airmass = np.atleast_1d(np.asarray(airmass, dtype=float))
    _color = np.atleast_1d(np.asarray(color, dtype=float))
    values = {
        'zp': np.ones_like(_airmass), 'airmass': _airmass, 'color': _color,
        'color2': _color**2, 'airmass_color': _airmass * _color
    }
    return np.stack([values[t] for t in terms])


@dataclass
class CalibrationFit:
    """
    The result of fitting a calibration model; see :func:`fit_calibration`.
    """
    terms: list
    """The names of the terms."""
    coeff: np.ndarray
    """The coefficients of all of the terms, including any fixed terms."""
    cov: np.ndarray
    """The covariance of the coefficients; zero for the fixed terms."""
    free: np.ndarray
    """Boolean array selecting the fitted (not fixed) terms."""
    gpm: np.ndarray
    """Boolean array selecting the measurements used in the fit."""
    chi2: float
    """The chi-square of the fit."""
    dof: int
    """The number of degrees of freedom of the fit."""
    rms: float
    """The root-mean-square residual of the fit."""

    @property
    def rchi2(self):
        """
        The reduced chi-square, or None if the fit has no degrees of freedom.
        """
        return self.chi2 / self.dof if self.dof > 0 else None

    @property
    def errors(self):
        """
        The errors of the coefficients; zero for the fixed terms.
        """
        return np.sqrt(np.diag(self.cov))

    def model(self, airmass, color):
        """
        Return the model difference between the instrumental and calibrated
        magnitudes, :math:`m - M`.
        """
        return np.dot(self.coeff, design_matrix(self.terms, airmass, color))

    def model_error(self, airmass, color):
        """
        Return the error of :meth:`model`, propagated from the covariance of
        the coefficients.
        """
        x = design_matrix(self.terms, airmass, color)
        return np.sqrt(np.sum(x * np.dot(self.cov, x), axis=0))


def fit_calibration(mag_c, color, mag_i, airmass, terms, err=None, fixed=None, gpm=None):
    r"""
    Fit a calibration model.

    The approach is the weighted linear least-squares method of Harris et al.
    (1981).

    Parameters
    ----------
    mag_c : :class:`numpy.ndarray`
        The catalog magnitudes, :math:`M`.
    color : :class:`numpy.ndarray`
        The catalog colors, :math:`C`.
    mag_i : :class:`numpy.ndarray`
        The instrumental magnitudes, :math:`m`.
    airmass : :class:`numpy.ndarray`
        The airmass of each measurement.
    terms : list
        The names of the terms; see :func:`check_terms`.
    err : :class:`numpy.ndarray`, optional
        The error of each measurement of :math:`m - M`, used to weight the
        fit.  If None, the weights are uniform.
    fixed : dict, optional
        The terms to fix, mapped to their values.
    gpm : :class:`numpy.ndarray`, optional
        Boolean array selecting the measurements to use.  If None, all are
        used.

    Returns
    -------
    :class:`CalibrationFit`
        The result of the fit.  If the measurements have errors and the
        reduced chi-square is larger than 1, the covariance is scaled by it.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if a fixed term is not in the model, there are no more
        measurements than fitted terms, or the fitted terms are degenerate.
    """
    _terms = check_terms(terms)
    _fixed = {} if fixed is None else fixed
    unknown = [t for t in _fixed if t not in _terms]
    if len(unknown) > 0:
        raise PipnickParameterError(f'Fixed terms are not in the model: {", ".join(unknown)}.')
    _mag_i = np.asarray(mag_i, dtype=float)
    _gpm = np.ones(_mag_i.shape, dtype=bool) if gpm is None else np.asarray(gpm, dtype=bool)
    _err = np.ones(_mag_i.shape, dtype=float) if err is None else np.asarray(err, dtype=float)
    free = np.array([t not in _fixed for t in _terms])
    nfree = int(np.sum(free))
    ngood = int(np.sum(_gpm))
    if ngood <= nfree:
        raise PipnickParameterError(
            f'Cannot fit {nfree} terms to {ngood} measurements; more measurements than terms '
            'are required.'
        )

    x = design_matrix(_terms, np.asarray(airmass)[_gpm], np.asarray(color)[_gpm])
    coeff = np.array([_fixed.get(t, 0.) for t in _terms], dtype=float)
    # The measurements, after subtracting the fixed terms
    y = _mag_i[_gpm] - np.asarray(mag_c, dtype=float)[_gpm] - np.dot(coeff, x)
    wgt = 1. / _err[_gpm]**2
    cov = np.zeros((len(_terms), len(_terms)), dtype=float)
    if nfree > 0:
        xf = x[free]
        # The fitted terms must be independent, e.g., the zero point and the
        # extinction cannot both be fit to measurements at a single airmass
        if np.linalg.matrix_rank(xf.T * np.sqrt(wgt)[:,None]) < nfree:
            raise PipnickParameterError(
                'The fitted terms are degenerate for these measurements (e.g., the airmass or '
                'color does not vary); fix or remove terms.'
            )
        # Eqns. 3.8-3.10 of Harris et al. (1981)
        coeff[free], cov[np.ix_(free, free)], _, _ = weighted_lstsq(
            xf.T, y, err=_err[_gpm], scale_cov=err is not None
        )
    resid = y - np.dot(coeff[free], x[free]) if nfree > 0 else y
    chi2 = float(np.sum(wgt * resid**2))
    dof = ngood - nfree
    return CalibrationFit(
        terms=_terms, coeff=coeff, cov=cov, free=free, gpm=_gpm, chi2=chi2, dof=dof,
        rms=float(np.sqrt(np.mean(resid**2)))
    )


def reject_outliers(fit, mag_c, color, mag_i, airmass, err=None, rej=3.):
    """
    Identify the outliers of a calibration fit.

    The residuals of the measurements used in the fit, normalized by their
    errors, are sigma-clipped once about their median at ``rej`` times their
    standard deviation; see :func:`astropy.stats.sigma_clip`.

    Parameters
    ----------
    fit : :class:`CalibrationFit`
        The fit.
    mag_c, color, mag_i, airmass, err : :class:`numpy.ndarray`
        The measurements; see :func:`fit_calibration`.
    rej : float, optional
        The rejection threshold, in standard deviations.

    Returns
    -------
    :class:`numpy.ndarray`
        Boolean array selecting the measurements newly identified as
        outliers.
    """
    _err = np.ones(np.shape(mag_i), dtype=float) if err is None else np.asarray(err, dtype=float)
    resid = (
        np.asarray(mag_i, dtype=float) - np.asarray(mag_c, dtype=float)
        - fit.model(airmass, color)
    ) / _err
    outlier = np.zeros(np.shape(mag_i), dtype=bool)
    if np.sum(fit.gpm) < 3:
        return outlier
    clipped = stats.sigma_clip(
        resid[fit.gpm], sigma=rej, maxiters=1, cenfunc='median', stdfunc='std'
    )
    outlier[np.where(fit.gpm)[0][np.ma.getmaskarray(clipped)]] = True
    return outlier


def zero_point(terms, coeff, cov, airmass):
    r"""
    Compute the zero point and color coefficients of an image.

    For an image at airmass :math:`X`, the calibrated magnitude of a source
    with color :math:`C` and count rate :math:`r` (electrons per second) is

    .. math::

        M = -2.5\log(r) + Z + c_1\ C + c_2\ C^2,

    with :math:`Z = -(k_{\rm zp} + k_{\rm airmass} X)`, :math:`c_1 =
    -(k_{\rm color} + k_{\rm airmass\_color} X)`, and :math:`c_2 = -k_{\rm
    color2}`, where any term not in the model is 0.

    Parameters
    ----------
    terms : list
        The names of the terms.
    coeff : :class:`numpy.ndarray`
        The coefficients.
    cov : :class:`numpy.ndarray`
        The covariance of the coefficients.
    airmass : float
        The airmass of the image.

    Returns
    -------
    dict
        The zero point (``zp``) and its statistical error (``zp_err``), and
        the color coefficients (``c1`` and ``c2``).
    """
    k = dict(zip(terms, np.asarray(coeff, dtype=float)))
    grad = np.array([-1. if t == 'zp' else (-airmass if t == 'airmass' else 0.) for t in terms])
    return {
        'zp': -(k['zp'] + k.get('airmass', 0.) * airmass),
        'zp_err': float(np.sqrt(np.dot(grad, np.dot(np.asarray(cov, dtype=float), grad)))),
        # NOTE: Subtracting from 0 avoids writing -0 for missing terms
        'c1': 0. - (k.get('color', 0.) + k.get('airmass_color', 0.) * airmass),
        'c2': 0. - k.get('color2', 0.)
    }


def residual_drift(time, resid, err):
    """
    Measure the linear trend of calibration residuals with time.

    This assesses a drift of the zero point through the night, which the
    calibration model does not include.

    Parameters
    ----------
    time : :class:`numpy.ndarray`
        The time of each measurement, e.g., in hours.
    resid : :class:`numpy.ndarray`
        The residuals (magnitudes).
    err : :class:`numpy.ndarray`
        The errors of the residuals, used to weight the fit.

    Returns
    -------
    tuple
        The slope (magnitudes per unit time) and its error, which is scaled by
        the square root of the reduced chi-square if that is larger than 1, or
        None if there are fewer than three measurements or they are all at
        the same time.
    """
    _time = np.asarray(time, dtype=float)
    if len(_time) < 3 or np.ptp(_time) == 0:
        return None
    x = np.column_stack([np.ones_like(_time), _time - np.mean(_time)])
    coeff, cov, _, _ = weighted_lstsq(x, resid, err=err)
    return float(coeff[1]), float(np.sqrt(cov[1,1]))


def iterate_colors(
    mag_i, err_i, airmass, source, band, calibrations, default_color=None, atol=1e-3,
    maxiter=20
):
    """
    Calibrate instrumental magnitudes, iterating to find the colors of the
    sources.

    The calibration of each band depends on the color of the source, which is
    computed from its calibrated magnitudes.  The colors start at their
    default values; each iteration calibrates all of the measurements with
    the current colors, averages the calibrated magnitudes of each source in
    each band (weighted by the instrumental errors, so that repeated images
    of a source in a band are combined), and recomputes the colors.  The
    iterations stop when the colors change by less than ``atol``.

    A color is computed only if the source has calibrated measurements in
    both of its bands; otherwise, the default color is used and flagged as
    assumed.

    Parameters
    ----------
    mag_i : :class:`numpy.ma.MaskedArray`
        The instrumental magnitude of each measurement; masked measurements
        are not calibrated.
    err_i : :class:`numpy.ma.MaskedArray`
        The errors of the instrumental magnitudes.  Masked errors are taken to
        be 1 when averaging.
    airmass : :class:`numpy.ndarray`
        The airmass of each measurement.
    source : :class:`numpy.ndarray`
        An identifier of the source of each measurement.
    band : :class:`numpy.ndarray`
        The band of each measurement.
    calibrations : dict
        The calibration of each band: a dictionary with the ``terms``, the
        coefficients (``coeff``), their covariance (``cov``), and the two
        bands of the color (``color``, e.g., ``('B', 'V')``).  Measurements in
        bands without a calibration are not calibrated.
    default_color : dict, optional
        The default value of each color, keyed by its bands, e.g.,
        ``{('B', 'V'): 0.6}``.  Colors not provided default to 0.
    atol : float, optional
        The convergence tolerance for the colors (magnitudes).
    maxiter : int, optional
        The maximum number of iterations.

    Returns
    -------
    mag, mag_err : :class:`numpy.ma.MaskedArray`
        The calibrated magnitudes and their statistical errors (the
        instrumental error and the error from the covariance of the
        coefficients), masked for measurements that are not calibrated (and,
        for the errors, those without an instrumental error).
    color : :class:`numpy.ma.MaskedArray`
        The color used to calibrate each measurement.
    assumed : :class:`numpy.ndarray`
        Boolean array flagging measurements calibrated with the default
        color.
    """
    _mag = np.ma.asarray(mag_i, dtype=float)
    _err = np.ma.filled(np.ma.asarray(err_i, dtype=float), 1.)
    _airmass = np.asarray(airmass, dtype=float)
    _source = np.asarray(source)
    _band = np.asarray(band)
    _default = {} if default_color is None else default_color
    n = len(_mag)
    valid = np.logical_not(np.ma.getmaskarray(_mag)) & np.isin(_band, list(calibrations.keys()))
    mag = np.zeros(n, dtype=float)
    mag_err = np.zeros(n, dtype=float)
    used_color = np.zeros(n, dtype=float)
    assumed = np.ones(n, dtype=bool)

    # The current color of each source, for each color used
    colors = {}
    for i in np.where(valid)[0]:
        pair = tuple(calibrations[_band[i]]['color'])
        colors[(_source[i], pair)] = _default.get(pair, 0.)
    known = {key: False for key in colors}

    for _ in range(maxiter):
        for i in np.where(valid)[0]:
            cal = calibrations[_band[i]]
            used_color[i] = colors[(_source[i], tuple(cal['color']))]
            x = design_matrix(cal['terms'], _airmass[i], used_color[i])[:,0]
            mag[i] = _mag[i] - np.dot(cal['coeff'], x)
        # Weighted mean magnitude of each source in each band
        means = {}
        for key in set(zip(_source[valid], _band[valid])):
            indx = valid & (_source == key[0]) & (_band == key[1])
            wgt = 1. / _err[indx]**2
            means[key] = np.sum(wgt * mag[indx]) / np.sum(wgt)
        change = 0.
        for key in colors:
            src, pair = key
            if (src, pair[0]) in means and (src, pair[1]) in means:
                new = means[(src, pair[0])] - means[(src, pair[1])]
                change = max(change, abs(new - colors[key]))
                colors[key] = new
                known[key] = True
        if change < atol:
            break

    for i in np.where(valid)[0]:
        cal = calibrations[_band[i]]
        key = (_source[i], tuple(cal['color']))
        assumed[i] = not known[key]
        used_color[i] = colors[key]
        x = design_matrix(cal['terms'], _airmass[i], used_color[i])
        mag[i] = _mag[i] - np.dot(cal['coeff'], x[:,0])
        mag_err[i] = np.sqrt(_err[i]**2 + np.dot(x[:,0], np.dot(cal['cov'], x[:,0])))
    invalid = np.logical_not(valid)
    no_err = invalid | np.ma.getmaskarray(np.ma.asarray(err_i))
    return (
        np.ma.MaskedArray(mag, mask=invalid), np.ma.MaskedArray(mag_err, mask=no_err),
        np.ma.MaskedArray(used_color, mask=invalid), assumed & valid
    )
