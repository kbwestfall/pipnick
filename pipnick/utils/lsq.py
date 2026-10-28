"""
Weighted linear least-squares fitting.
"""

import numpy as np

from scipy import linalg


def weighted_lstsq(x, y, err=None, gpm=None, scale_cov=True):
    r"""
    Solve a weighted linear least-squares problem by its normal equations.

    The model is :math:`y = X c`, and the solution minimizes
    :math:`\chi^2 = \sum_k w_k (y_k - \sum_i X_{ki} c_i)^2`, with weights
    :math:`w_k = 1/\sigma_k^2`.  The normal equations,
    :math:`(X^T W X) c = X^T W y`, are solved directly, and the covariance of
    the coefficients is :math:`(X^T W X)^{-1}`.

    Parameters
    ----------
    x : :class:`numpy.ndarray`
        The design matrix, :math:`X`, with shape ``(nmeas, npar)``.
    y : :class:`numpy.ndarray`
        The measurements, with shape ``(nmeas,)``.
    err : :class:`numpy.ndarray`, optional
        The error of each measurement, :math:`\sigma_k`.  If None, the
        weights are uniform.
    gpm : :class:`numpy.ndarray`, optional
        Boolean array selecting the measurements to use.  If None, all are
        used.
    scale_cov : bool, optional
        Scale the covariance by the reduced :math:`\chi^2`, if that is larger
        than 1 (and there are more measurements than parameters).

    Returns
    -------
    coeff : :class:`numpy.ndarray`
        The coefficients, :math:`c`.
    cov : :class:`numpy.ndarray`
        The covariance of the coefficients.
    chi2 : float
        The :math:`\chi^2` of the fit, over the selected measurements.
    dof : int
        The number of degrees of freedom: the number of selected measurements
        minus the number of parameters.

    Raises
    ------
    :exc:`numpy.linalg.LinAlgError`
        Raised if the normal equations are singular.
    """
    _x = np.asarray(x, dtype=float)
    _y = np.asarray(y, dtype=float)
    _gpm = np.ones(_y.shape, dtype=bool) if gpm is None else np.asarray(gpm, dtype=bool)
    _x = _x[_gpm]
    _y = _y[_gpm]
    if err is None:
        wgt = np.ones(_y.shape, dtype=float)
    else:
        wgt = 1. / np.asarray(err, dtype=float)[_gpm]**2
    normal = np.dot(_x.T * wgt, _x)
    coeff = linalg.solve(normal, np.dot(_x.T * wgt, _y))
    cov = linalg.inv(normal)
    chi2 = float(np.sum(wgt * (_y - np.dot(_x, coeff))**2))
    dof = len(_y) - _x.shape[1]
    if scale_cov and dof > 0 and chi2 / dof > 1.:
        cov *= chi2 / dof
    return coeff, cov, chi2, dof
