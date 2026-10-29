"""
Tests for the detection of sources and the background of an image.
"""

import numpy as np
import pytest

from pipnick.astrometry.cache import anet_input_hash
from pipnick.image.background import estimate_background, global_background, masked_median
from pipnick.image.detect import SourceFinder, brightest, find_sources, isolated
from pipnick.pipelines.astrometry import _detect
from pipnick.pkg.exceptions import PipnickParameterError

SHAPE = (200, 220)
FWHM = 4.


def star_image(seed=3, level=100., noise=3.):
    """
    Render an image of 30 Gaussian stars on a constant background, with a
    masked column of bad pixels.
    """
    rng = np.random.default_rng(seed)
    x = rng.uniform(10., SHAPE[1] - 10., 30)
    y = rng.uniform(10., SHAPE[0] - 10., 30)
    flux = 10**rng.uniform(3., 4.5, 30)
    yy, xx = np.mgrid[:SHAPE[0], :SHAPE[1]]
    sig = FWHM / 2.355
    data = rng.normal(level, noise, size=SHAPE)
    for xi, yi, fi in zip(x, y, flux):
        data += fi / (2 * np.pi * sig**2) * np.exp(-((xx - xi)**2 + (yy - yi)**2) / (2 * sig**2))
    mask = np.zeros(SHAPE, dtype=bool)
    mask[:,100] = True
    data[mask] = 1e5
    return data, mask, (x, y, flux)


def test_find_sources():
    data, mask, (x, y, flux) = star_image()
    sources = find_sources(data, mask=mask, fwhm=FWHM, thresh=10.)
    assert sources.meta['background'] == pytest.approx(100., abs=0.5), (
        'The background level should be measured from the unmasked pixels'
    )
    assert sources.meta['noise'] == pytest.approx(3., rel=0.1), (
        'The noise should be measured from the unmasked pixels'
    )
    dist = np.hypot(x[:,None] - sources['x'][None,:], y[:,None] - sources['y'][None,:])
    assert np.sum(np.min(dist, axis=0) < 1.) >= 0.9 * len(sources), (
        'Nearly all detections should be the injected stars'
    )
    assert np.all(np.absolute(np.asarray(sources['x']) - 100.) > 1.), (
        'The masked column should not be detected'
    )

    # The masked pixels are ignored, whatever their values
    masked = np.ma.MaskedArray(np.where(mask, -5., data), mask=mask)
    other = find_sources(masked, fwhm=FWHM, thresh=10.)
    assert np.array_equal(other['x'], sources['x']), (
        'The values of the masked pixels should not matter with masked=background'
    )
    zero = find_sources(data, mask=mask, fwhm=FWHM, thresh=10., masked='zero')
    assert len(zero) > 0, 'Setting the masked pixels to 0 should still detect the stars'

    empty = find_sources(np.full(SHAPE, 100.), fwhm=FWHM)
    assert len(empty) == 0 and empty.colnames[:3] == ['x', 'y', 'flux'], (
        'An image without noise should give an empty table with the same columns'
    )
    with pytest.raises(PipnickParameterError):
        find_sources(data, background='none')
    with pytest.raises(PipnickParameterError):
        find_sources(data, masked='nan')


def test_source_finder():
    data, mask, _ = star_image()
    sources = find_sources(data, mask=mask, fwhm=FWHM, thresh=10.)
    finder = SourceFinder(FWHM, 10. * sources.meta['noise'])
    found = finder(np.where(mask, 0., data - sources.meta['background']), mask=mask)
    assert np.array_equal(found['x'], sources['x']), (
        'The finder should detect the same sources at the same threshold'
    )
    assert finder(np.zeros(SHAPE), mask=mask) is None, (
        'The finder should return None if no sources are found, as photutils expects'
    )


def test_selections():
    data, mask, _ = star_image()
    sources = find_sources(data, mask=mask, fwhm=FWHM, thresh=10.)
    top = brightest(sources, 5)
    assert len(top) == 5 and np.all(np.diff(top['flux']) <= 0), (
        'The brightest sources should be ordered by decreasing flux'
    )
    assert top['flux'][0] == np.max(sources['flux']), 'The brightest source should be first'
    assert len(brightest(sources)) == len(sources), 'All sources should be kept without n'

    x = np.array([0., 3., 50.])
    y = np.zeros(3)
    flux = np.ma.MaskedArray([10., 20., 5.])
    assert list(isolated(x, y, flux, 5.)) == [False, True, True], (
        'Only the brightest of neighbors, and sources without neighbors, should be selected'
    )


def test_astrometry_detection():
    data, mask, _ = star_image()
    data[mask] = 0.
    sources = _detect(data, mask, FWHM, 5., 20)
    assert sources.colnames == ['x', 'y', 'flux'] and len(sources) == 20, (
        'The astrometry should detect the brightest sources'
    )
    # The astrometry.net cache key of the sources of this image, as computed by
    # the detection of the astrometry before it used find_sources.  A change
    # would mean that the astrometry's detections, or their ranking by the
    # flux of DAOStarFinder, changed (as they did in photutils 3.1, which
    # changed the flux to the sum within the kernel footprint), so that cached
    # astrometry.net results would no longer be found and images could be
    # submitted again.
    hints = {'scale_units': 'arcsecperpix', 'scale_lower': 0.3, 'scale_upper': 0.4}
    assert anet_input_hash(data, sources, hints, False) == '7117b0f8ff89602e', (
        'The detections of the astrometry, and so the astrometry.net cache keys, should '
        'not change'
    )


def test_backgrounds():
    data, mask, _ = star_image()
    level, noise = global_background(data, mask=mask)
    assert level == pytest.approx(100., abs=0.5) and noise == pytest.approx(3., rel=0.1), (
        'The global background should be measured from the unmasked pixels'
    )
    assert masked_median(data, mask=mask) == pytest.approx(100., abs=0.2), (
        'The median with the sources masked should be the background level'
    )
    assert estimate_background(data, mask=mask, method='none') == 0., (
        'No background should be 0'
    )
    bkg = estimate_background(data, mask=mask, method='2d', box_size=50)
    assert bkg.shape == SHAPE and np.allclose(np.median(bkg), 100., atol=0.5), (
        'The 2D background should be the background level'
    )
    with pytest.raises(PipnickParameterError):
        estimate_background(data, method='bad')
