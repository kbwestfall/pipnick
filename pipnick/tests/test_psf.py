"""
Tests for the PSF model.
"""

import numpy as np
import pytest

from pipnick.datamodel.psf import PSFModel
from pipnick.photometry.moffat import FitMoffat2D, fwhm_to_gamma, gamma_to_fwhm
from pipnick.pkg.exceptions import PipnickParameterError

# An elliptical model: x0, y0, amplitude, gamma1, gamma2, phi, alpha, background
ELLIP = [10., 12., 100., 2.0, 3.0, 0.4, 2.8, 1.]


def test_conversions():
    assert np.isclose(gamma_to_fwhm(fwhm_to_gamma(4.2, 2.5), 2.5), 4.2), (
        'The FWHM and core-width conversions should be inverses'
    )
    assert np.isclose(gamma_to_fwhm(1., 1e6), 0., atol=1e-2), (
        'A Moffat function with a large power index should be narrow for its core width'
    )


def test_fwhm_definitions():
    psf = PSFModel('ellip', ELLIP)
    fwhm = [2 * g * np.sqrt(2**(1 / ELLIP[6]) - 1) for g in ELLIP[3:5]]
    assert np.isclose(psf.fwhm, np.mean(fwhm), rtol=1e-14, atol=0.), (
        'The FWHM should be the mean of the FWHM along the two axes (the definition of the '
        'photometry tables)'
    )
    assert np.isclose(psf.fwhm_major, max(fwhm), rtol=1e-14, atol=0.), (
        'The major-axis FWHM should be the larger of the two (used by the PSF matching)'
    )
    assert np.isclose(psf.fwhm_minor, min(fwhm), rtol=1e-14, atol=0.), (
        'The minor-axis FWHM should be the smaller of the two'
    )
    assert psf.beta == ELLIP[6] and psf.phi == ELLIP[5], (
        'The power index and rotation should come from the parameters'
    )
    circ = PSFModel.circular(4.5, 3.)
    assert np.isclose(circ.fwhm, 4.5) and circ.fwhm == circ.fwhm_major == circ.fwhm_minor, (
        'All of the FWHM definitions should agree for a circular model'
    )


def test_construction():
    with pytest.raises(PipnickParameterError):
        PSFModel('gauss', ELLIP)
    with pytest.raises(PipnickParameterError):
        PSFModel('circ', ELLIP)
    fit = FitMoffat2D(PSFModel.circular(4., 2.5).render(21))
    assert PSFModel.from_fit(fit).fittype == 'circ', (
        'A model from a circular fit should be circular'
    )
    psf = PSFModel('ellip', ELLIP)
    meta = psf.to_meta(0.333, 25)
    assert meta['fwhm_pix'] == psf.fwhm and np.isclose(meta['fwhm_arcsec'], psf.fwhm * 0.333), (
        'The metadata should record the FWHM in pixels and arcsec'
    )
    assert PSFModel.from_meta(meta).par == psf.par, (
        'A model should be recovered from its metadata'
    )


def test_render():
    for psf in [PSFModel.circular(4., 2.5), PSFModel('ellip', ELLIP)]:
        image = psf.render(31)
        assert np.isclose(np.sum(image), 1.), 'The rendered PSF should be normalized'
        assert np.unravel_index(np.argmax(image), image.shape) == (15, 15), (
            'The rendered PSF should be centered'
        )
    # The photutils model has the same profile
    psf = PSFModel.circular(4., 2.5)
    model = psf.to_photutils(flux=1., x_0=15., y_0=15.)
    yy, xx = np.mgrid[:31, :31]
    expected = model(xx, yy)
    assert np.allclose(psf.render(31), expected / np.sum(expected)), (
        'The photutils model should have the shape of the PSF model'
    )
    # The number of pixels above half the peak along the central row
    narrow = psf.render(31)
    wide = psf.render(31, scale=2.)
    width = np.sum(narrow[15] > narrow[15,15] / 2)
    assert abs(np.sum(wide[15] > wide[15,15] / 2) - 2 * width) <= 1, (
        'Rendering at twice the scale should double the width of the PSF'
    )
