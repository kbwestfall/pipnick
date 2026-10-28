"""
Tests for the weighted linear least-squares fit.
"""

import numpy as np

from pipnick.utils.lsq import weighted_lstsq


def line(n, seed=1):
    """A straight line with Gaussian errors of varying size."""
    rng = np.random.default_rng(seed)
    t = np.linspace(-1., 1., n)
    err = rng.uniform(0.5, 2., n)
    y = 3. + 2. * t + rng.normal(0., err)
    return np.column_stack([np.ones(n), t]), y, err


def test_matches_numpy():
    x, y, err = line(50)
    coeff, cov, chi2, dof = weighted_lstsq(x, y, err=err, scale_cov=False)
    expected = np.linalg.lstsq(x / err[:,None], y / err, rcond=None)[0]
    assert np.allclose(coeff, expected, rtol=1e-12, atol=0.), (
        'The solution should match the weighted solution of numpy.linalg.lstsq'
    )
    expected = np.linalg.inv(np.dot(x.T / err**2, x))
    assert np.allclose(cov, expected, rtol=1e-12, atol=0.), (
        'The covariance should be the inverse of the weighted normal matrix'
    )
    assert dof == 48, 'The degrees of freedom should be the measurements minus the parameters'
    assert np.isclose(chi2, np.sum(((y - np.dot(x, coeff)) / err)**2)), (
        'The chi-square should be the weighted sum of the squared residuals'
    )


def test_selection_and_scaling():
    x, y, err = line(50)
    gpm = np.ones(50, dtype=bool)
    gpm[::7] = False
    y[::7] += 100.
    coeff, cov, chi2, dof = weighted_lstsq(x, y, err=err, gpm=gpm, scale_cov=False)
    _coeff, _cov, _, _ = weighted_lstsq(x[gpm], y[gpm], err=err[gpm], scale_cov=False)
    assert np.array_equal(coeff, _coeff) and np.array_equal(cov, _cov), (
        'Deselected measurements should be ignored'
    )
    # Errors that are too small by a factor of 3 inflate the chi-square by 9
    coeff, cov, chi2, dof = weighted_lstsq(x, y, err=err / 3., gpm=gpm, scale_cov=False)
    _coeff, scaled, _, _ = weighted_lstsq(x, y, err=err / 3., gpm=gpm)
    assert np.allclose(scaled, cov * chi2 / dof), (
        'The covariance should be scaled by the reduced chi-square when it exceeds 1'
    )
    assert chi2 / dof > 1., 'The underestimated errors should give a reduced chi-square above 1'


def test_uniform_weights():
    x, y, _ = line(20)
    coeff, _, _, _ = weighted_lstsq(x, y)
    expected = np.linalg.lstsq(x, y, rcond=None)[0]
    assert np.allclose(coeff, expected, rtol=1e-12, atol=0.), (
        'Without errors, the fit should be the ordinary least-squares solution'
    )
