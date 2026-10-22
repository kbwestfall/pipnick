"""
Tests for the camera definitions, using synthetic headers.
"""

import numpy as np
import pytest
from astropy.io import fits

from pipnick import cameras
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickParameterError
from pipnick.utils import masking


def c2_header(**kwargs):
    """
    Construct a minimal header for the original Nickel camera.
    """
    hdr = fits.Header({'NAXIS1': 1056, 'NAXIS2': 1024, 'INSTRUME': 'Nickel Direct Camera',
                       'DSENSOR': 'Loral 2Kx2K', 'OBSTYPE': 'OBJECT', 'OBJECT': 'NGC_6543',
                       'EXPTIME': 60, 'FILTNAM': 'R', 'CDELT1U': -2, 'CDELT2U': -2,
                       'READ-SPD': 80, 'GAIN': 1, 'DATASEC': '[1:1024,1:1024]', 'COVER': 32})
    hdr.update(kwargs)
    return hdr


def andor_header(binning=2, window=None, **kwargs):
    """
    Construct a minimal header for the Nickel Andor camera.
    """
    x1, x2, y1, y2 = (1, 2088, 1, 2048) if window is None else window
    hdr = fits.Header({'NAXIS1': (x2 - x1 + 1) // binning, 'NAXIS2': (y2 - y1 + 1) // binning,
                       'CAMID': 'headmodel DZ936_BV serial 31382', 'INSTRUME': 'Direct',
                       'OBSTYPE': '?', 'OBJECT': 'M11 V', 'XPOSURE': 30.005, 'FILTER': 'V',
                       'CDELT1U': binning, 'CDELT2U': binning, 'AMPCONF': '1.0MHz',
                       'WINDOW': f'{x1} {x2} {y1} {y2}'})
    hdr.update(kwargs)
    return hdr


def test_registry():
    assert sorted(cameras.available_cameras()) == ['NickelAndor', 'NickelC2']
    assert cameras.get_camera('NickelC2').name == 'NickelC2'
    with pytest.raises(PipnickParameterError):
        cameras.get_camera('NotACamera')


def test_bad_config():
    with pytest.raises(PipnickParameterError):
        cameras.NickelCamera({'camera': {'name': 'Test'}})


def test_identify():
    assert cameras.identify_camera(c2_header()).name == 'NickelC2'
    assert cameras.identify_camera(andor_header()).name == 'NickelAndor'
    # Neither camera
    with pytest.raises(PipnickDataModelError):
        cameras.identify_camera(fits.Header({'INSTRUME': 'Direct'}))
    # Both cameras
    with pytest.raises(PipnickDataModelError):
        cameras.identify_camera(c2_header(CAMID='headmodel DZ936_BV serial 31382'))


def test_keywords():
    c2 = cameras.get_camera('NickelC2')
    andor = cameras.get_camera('NickelAndor')
    assert c2.header_value(c2_header(), 'filter') == 'R'
    assert c2.header_value(c2_header(), 'exptime') == 60
    assert andor.header_value(andor_header(), 'filter') == 'V'
    assert andor.header_value(andor_header(), 'exptime') == 30.005
    assert andor.header_value(andor_header(), 'datasec') is None
    # Standard keywords take precedence when present
    hdr = andor_header(EXPTIME=30.0, FILTNAM='R')
    assert andor.header_value(hdr, 'exptime') == 30.0
    assert andor.header_value(hdr, 'filter') == 'R'


def test_binning():
    c2 = cameras.get_camera('NickelC2')
    andor = cameras.get_camera('NickelAndor')
    # The old camera reports negative binning
    assert c2.binning(c2_header()) == (2, 2)
    assert andor.binning(andor_header(binning=4)) == (4, 4)
    # CCDSUM takes precedence
    assert andor.binning(andor_header(CCDSUM='1 1')) == (1, 1)
    with pytest.raises(PipnickDataModelError):
        andor.binning(fits.Header({'CAMID': 'headmodel DZ936'}))


def test_readout():
    c2 = cameras.get_camera('NickelC2')
    andor = cameras.get_camera('NickelAndor')

    readout = c2.readout(c2_header())
    assert readout.mode == 'slow' and readout.gain == 1.8 and readout.readnoise == 10.7
    assert not readout.placeholder
    # The old camera's GAIN keyword is a DCS gain index; it must not be used
    assert c2.readout(c2_header(GAIN=1)).gain == 1.8
    # Undefined readout settings
    assert c2.readout(c2_header(**{'READ-SPD': 10})).mode == 'unknown'

    assert andor.readout(andor_header()).mode == 'fast'
    readout = andor.readout(andor_header(AMPCONF='0.05MHz'))
    assert readout.mode == 'slow' and readout.placeholder
    assert readout.saturation == andor.saturation
    # Header values take precedence and are not placeholders
    readout = andor.readout(andor_header(GAIN=2.5, RDNOISE=4.0))
    assert readout.gain == 2.5 and readout.readnoise == 4.0 and not readout.placeholder


def test_c2_sections():
    c2 = cameras.get_camera('NickelC2')
    hdr = c2_header()
    # From the header
    assert c2.sections(hdr) == ('[1:1024,1:1024]', '[1025:1056,1:1024]')
    # From the configured geometry
    del hdr['DATASEC']
    assert c2.sections(hdr) == ('[1:1024,1:1024]', '[1025:1056,1:1024]')


@pytest.mark.parametrize('binning, datasec, biassec', [
    (1, '[20:2067,1:2048]', '[2075:2088,1:2048]'),
    (2, '[11:1033,1:1024]', '[1038:1044,1:1024]'),
    (4, '[6:516,1:512]', '[520:522,1:512]'),
])
def test_andor_sections(binning, datasec, biassec):
    # These sections were measured from the 2026 commissioning data
    andor = cameras.get_camera('NickelAndor')
    assert andor.sections(andor_header(binning=binning)) == (datasec, biassec)


def test_andor_window():
    andor = cameras.get_camera('NickelAndor')
    # The windowed frame from the 2026 commissioning run
    hdr = andor_header(window=(960, 1205, 840, 1069))
    assert (hdr['NAXIS1'], hdr['NAXIS2']) == (123, 115)
    assert andor.is_windowed(hdr)
    assert not andor.is_windowed(andor_header())
    # The window excludes the overscan
    assert andor.sections(hdr) == ('[1:123,1:115]', None)
    assert andor.bpm(hdr).shape == (115, 123)


def test_shape_mismatch():
    andor = cameras.get_camera('NickelAndor')
    with pytest.raises(PipnickDataModelError):
        andor.sections(andor_header(NAXIS1=1056))


def test_header_sections():
    # Sections in the header take precedence over the configured geometry
    andor = cameras.get_camera('NickelAndor')
    hdr = andor_header(DATASEC='[12:1032,1:1024]', BIASSEC='[1036:1044,1:1024]')
    assert andor.sections(hdr) == ('[12:1032,1:1024]', '[1036:1044,1:1024]')


def test_c2_bpm():
    c2 = cameras.get_camera('NickelC2')
    bpm = c2.bpm(c2_header())
    assert bpm.shape == (1024, 1056)
    # The legacy 2x2-binned mask, plus the always-saturated first and last data
    # columns and the half-illuminated first row
    legacy = masking.create_mask((1024, 1056), bad_cols=[0, 255, 256, 783, 784, 1002, 1023],
                                 bad_regions=[np.array([[-0.5, 33.5], [34.5, -0.5], [-0.5, -0.5]]),
                                              np.array([[-0.5, 960.5], [64.5, 1024.5],
                                                        [-0.5, 1024.5]])])
    legacy[0] = True
    # Binning the unbinned mask flags any binned pixel with a bad unbinned
    # pixel, so it adds pixels along the diagonal edges of the blind corners
    assert np.all(bpm[legacy])
    extra = bpm & ~legacy
    assert np.all(np.sum(extra, axis=1) <= 1)
    assert np.sum(extra) < 50


def test_andor_bpm():
    andor = cameras.get_camera('NickelAndor')
    for binning in [1, 2, 4]:
        bpm = andor.bpm(andor_header(binning=binning))
        assert bpm.shape == (2048 // binning, 2088 // binning)
        assert not np.any(bpm)


@pytest.mark.parametrize('kwargs, frametype', [
    ({'OBJECT': 'Bias', 'EXPTIME': 0}, 'bias'),
    # Biases taken with OBSTYPE = 'OBJECT'
    ({'OBJECT': 'Bias', 'EXPTIME': 0, 'OBSTYPE': 'OBJECT'}, 'bias'),
    ({'OBJECT': 'Post-ievent recovery. Test bias.', 'EXPTIME': 0}, 'bias'),
    # Mislabelled: an exposed frame with OBJECT = 'Bias'
    ({'OBJECT': 'Bias', 'EXPTIME': 10}, 'unknown'),
    ({'OBJECT': 'Flat', 'EXPTIME': 12}, 'skyflat'),
    ({'OBJECT': 'Dome flat', 'EXPTIME': 30}, 'domeflat'),
    ({'OBJECT': 'Focusing14', 'EXPTIME': 10}, 'focus'),
    ({'OBJECT': 'Focus test B', 'EXPTIME': 10}, 'focus'),
    ({'OBJECT': 'Point', 'EXPTIME': 10}, 'pointing'),
    ({'OBJECT': 'test', 'EXPTIME': 10}, 'test'),
    ({'OBJECT': '109_199', 'EXPTIME': 30}, 'science'),
])
def test_c2_classify(kwargs, frametype):
    c2 = cameras.get_camera('NickelC2')
    assert c2.classify(c2_header(**kwargs)) == frametype


@pytest.mark.parametrize('kwargs, frametype', [
    ({'OBSTYPE': 'DARK', 'OBJECT': 'bias', 'XPOSURE': 0.207}, 'bias'),
    # A closed-shutter bias with a stale OBJECT
    ({'OBSTYPE': 'DARK', 'OBJECT': 'test', 'XPOSURE': 0.207}, 'bias'),
    ({'OBSTYPE': 'DARK', 'OBJECT': 'Dark', 'XPOSURE': 120.006}, 'dark'),
    # Open shutter, but labelled as a bias
    ({'OBJECT': 'bias', 'XPOSURE': 0.2}, 'unknown'),
    ({'OBJECT': 'Dome Flat V', 'XPOSURE': 3.0}, 'domeflat'),
    ({'OBJECT': 'Domeflat V', 'XPOSURE': 3.0}, 'domeflat'),
    ({'OBJECT': 'twilight flat slow readout', 'XPOSURE': 30.0}, 'skyflat'),
    ({'OBJECT': 'Focus 18 V Window', 'XPOSURE': 5.0}, 'focus'),
    ({'OBJECT': 'M11 (focus)', 'XPOSURE': 5.0}, 'focus'),
    ({'OBJECT': 'pointing18', 'XPOSURE': 5.0}, 'pointing'),
    ({'OBJECT': 'test open', 'XPOSURE': 1.0}, 'test'),
    ({'OBJECT': 'Northern Local Void', 'XPOSURE': 120.0}, 'science'),
    # An explicit frame-type keyword takes precedence
    ({'OBJECT': 'Northern Local Void', 'XPOSURE': 3.0, 'IMAGETYP': 'domeflat'}, 'domeflat'),
    # ... unless its value is not recognized
    ({'OBJECT': 'Northern Local Void', 'XPOSURE': 3.0, 'IMAGETYP': 'other'}, 'science'),
])
def test_andor_classify(kwargs, frametype):
    andor = cameras.get_camera('NickelAndor')
    assert andor.classify(andor_header(**kwargs)) == frametype
