"""
Tests for :mod:`pipnick.utils.combiner`.

These tests use only the module, ``numpy``, and ``astropy``, so that they can
be contributed to ``ccdproc`` with the module.
"""

import numpy as np
import pytest
from astropy import units
from astropy.nddata import CCDData, StdDevUncertainty, VarianceUncertainty

from pipnick.utils.combiner import Combiner


def images(sigma, n=4, shape=(5, 6), seed=1):
    """
    Return ``n`` images with constant uncertainties ``sigma`` (a scalar or a
    sequence of length ``n``).
    """
    rng = np.random.default_rng(seed)
    _sigma = np.broadcast_to(np.atleast_1d(sigma), (n,))
    return [
        CCDData(
            rng.normal(10., 1., shape), unit=units.electron,
            uncertainty=StdDevUncertainty(np.full(shape, s))
        ) for s in _sigma
    ]


def test_propagated_average():
    combined = Combiner(images(2.)).average_combine(propagate_uncertainty=True)
    assert np.allclose(combined.uncertainty.array, 1.), (
        'The uncertainty of the mean of 4 images with sigma 2 should be 1'
    )
    combined = Combiner(images([1., 2., 3., 4.])).average_combine(propagate_uncertainty=True)
    assert np.allclose(combined.uncertainty.array, np.sqrt(1 + 4 + 9 + 16) / 4), (
        'Different uncertainties should be added in quadrature'
    )
    default = Combiner(images(2.)).average_combine()
    assert not np.allclose(default.uncertainty.array, 1.), (
        'Without propagation, the uncertainty should be the scatter of the values'
    )


def test_propagated_average_masked():
    ccds = images(2.)
    combiner = Combiner(ccds)
    # Reject the first image at one pixel, and all images at another
    combiner.data_arr.mask[0,1,1] = True
    combiner.data_arr.mask[:,2,2] = True
    combined = combiner.average_combine(propagate_uncertainty=True)
    assert np.isclose(combined.uncertainty.array[1,1], 2. / np.sqrt(3)), (
        'Only the images that are not masked should contribute'
    )
    assert combined.uncertainty.array[2,2] == 0. and combined.mask[2,2], (
        'Where all images are masked, the uncertainty should be 0 and the pixel masked'
    )
    assert np.isclose(combined.uncertainty.array[0,0], 1.), 'Other pixels should be unaffected'


def test_propagated_average_scaled_weighted():
    combiner = Combiner(images(2.))
    combiner.scaling = [1., 2., 1., 2.]
    combined = combiner.average_combine(propagate_uncertainty=True)
    assert np.allclose(combined.uncertainty.array, np.sqrt(2 * 4 + 2 * 16) / 4), (
        'The uncertainties should be scaled with the images'
    )
    combiner = Combiner(images(2.))
    combiner.weights = np.array([1., 1., 2., 0.])
    combined = combiner.average_combine(propagate_uncertainty=True)
    assert np.allclose(combined.uncertainty.array, np.sqrt(4 + 4 + 16) / 4), (
        'The uncertainties should be weighted like the images'
    )


def test_uncertainty_types():
    ccds = images(2.)
    ccds[0].uncertainty = VarianceUncertainty(np.full(ccds[0].shape, 4.))
    combined = Combiner(ccds).average_combine(propagate_uncertainty=True)
    assert np.allclose(combined.uncertainty.array, 1.), (
        'Variance uncertainties should be converted to standard deviations'
    )


def test_missing_uncertainty():
    ccds = images(2.)
    ccds[0].uncertainty = None
    combiner = Combiner(ccds)
    assert combiner.uncertainty_arr is None, 'The uncertainties cannot all be kept'
    with pytest.raises(ValueError, match='Cannot propagate'):
        combiner.average_combine(propagate_uncertainty=True)
    assert combiner.average_combine().uncertainty is not None, (
        'The default combination should still work'
    )


def test_propagation_matches_formula():
    # The unweighted propagation (by astropy) should match the explicit formula
    rng = np.random.default_rng(2)
    sigma = rng.uniform(1., 5., (4, 5, 6))
    ccds = [
        CCDData(rng.normal(10., 1., (5, 6)), unit=units.electron, uncertainty=StdDevUncertainty(s))
        for s in sigma
    ]
    combiner = Combiner(ccds)
    combiner.scaling = [1., 2., 0.5, 1.5]
    combiner.data_arr.mask = rng.uniform(size=(4, 5, 6)) < 0.3
    combined = combiner.average_combine(propagate_uncertainty=True)
    scale = np.array([1., 2., 0.5, 1.5])[:,None,None]
    good = np.logical_not(combiner.data_arr.mask)
    n = np.sum(good, axis=0)
    expected = np.sqrt(np.sum((good * scale * sigma)**2, axis=0)) / np.where(n > 0, n, 1)
    assert np.allclose(combined.uncertainty.array, expected), (
        'The propagated uncertainty should match the formula for the scaled, masked average'
    )
