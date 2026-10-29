"""
Tests for the settings of the processing steps.
"""

import dataclasses

import pytest

from pipnick.cameras import get_camera
from pipnick.pkg.exceptions import PipnickParameterError
from pipnick.pkg.provenance import json_hash
from pipnick.pkg.settings import (
    AstrometrySettings, CalibrationSettings, MosaicSettings, PhotometrySettings
)


def test_photometry_settings():
    settings = PhotometrySettings(thresh=8, fwhm=1.5, group=0)
    # The settings dictionary of the photometry before the dataclass, whose
    # hash is recorded in every source table
    before = {
        'thresh': 8., 'fwhm': 1.5, 'fittype': 'circ', 'mode': 'all', 'stack': 'mean',
        'group': False, 'min_sep': 1., 'group_sep': 2., 'aper_radius': 2., 'psf_thresh': 15.
    }
    assert settings.as_dict() == before, 'The settings should be recorded as before'
    assert settings.hash() == json_hash(before) == 'f16a0463a4fddc0d', (
        'The hash of the settings (settings_hash) should not change'
    )
    assert isinstance(settings.thresh, float) and settings.group is False, (
        'The values should be normalized'
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        settings.thresh = 5.
    with pytest.raises(PipnickParameterError, match='fittype'):
        PhotometrySettings(fittype='gauss')
    with pytest.raises(PipnickParameterError, match='thresh'):
        PhotometrySettings(thresh=0.)
    with pytest.raises(PipnickParameterError, match='Resolve'):
        PhotometrySettings().hash()
    camera = get_camera('NickelAndor')
    assert PhotometrySettings().resolve(camera).fwhm == camera.fwhm_guess, (
        'The FWHM should default to the typical value for the camera'
    )
    assert PhotometrySettings(fwhm=2.).resolve(camera).fwhm == 2., (
        'A given FWHM should not be changed'
    )


def test_calibration_settings():
    settings = CalibrationSettings(terms=('color', 'zp'), colors={'R': 'R-I'}, maxiter=2.)
    assert settings.terms == ['zp', 'color'], 'The terms should be in their standard order'
    assert settings.colors['R'] == 'R-I' and settings.colors['V'] == 'B-V', (
        'The given colors should be combined with the defaults'
    )
    assert isinstance(settings.maxiter, int), 'The values should be normalized'
    before = {
        'terms': ['zp', 'airmass', 'color'],
        'colors': {'U': 'U-B', 'B': 'B-V', 'V': 'B-V', 'R': 'V-R', 'I': 'V-I'}, 'mag': 'psf',
        'match_radius': 2., 'max_neighbor_flux': 0.01, 'residual_radius': 1.5,
        'max_cat_err': 0.05, 'min_cat_obs': 1, 'sys_err': 0.01, 'rej': 3., 'maxiter': 5,
        'min_stars': 3, 'min_airmass_span': 0.2, 'min_color_span': 0.3, 'force': False,
        'fallback': True, 'drift_warn': 0.02
    }
    assert CalibrationSettings().hash() == json_hash(before) == 'd7e8b1b506a9f718', (
        'The hash of the settings, part of the input hash of each night, should not change'
    )
    changed = dataclasses.replace(settings, maxiter=0)
    assert changed.maxiter == 0 and changed.terms == settings.terms, (
        'A replaced setting should keep the others'
    )
    with pytest.raises(PipnickParameterError, match='mag'):
        CalibrationSettings(mag='kron')
    with pytest.raises(PipnickParameterError):
        CalibrationSettings(terms=['airmass'])
    with pytest.raises(PipnickParameterError):
        CalibrationSettings(colors={'V': 'B+V'})


def test_mosaic_settings():
    before = {
        'psf_match': 'mosaic', 'kernel': 'moffat', 'psf_margin': 1.05, 'scale': 'sources',
        'reference': None, 'match_radius': 1., 'min_common': 3, 'min_snr': 20.,
        'background': 'constant', 'match_background': False, 'reproject': 'adaptive',
        'pixel_scale': None, 'rej': 3., 'rej_frac': 0.1
    }
    assert MosaicSettings().hash() == json_hash(before) == 'c8e17799a8331df0', (
        'The hash of the settings, part of the hash of each mosaic, should not change'
    )
    assert MosaicSettings(rej=None).rej is None, 'No rejection should be allowed'
    with pytest.raises(PipnickParameterError, match='psf_match'):
        MosaicSettings(psf_match='image')
    with pytest.raises(PipnickParameterError, match='reproject'):
        MosaicSettings(reproject='nearest')


def test_astrometry_settings():
    settings = AstrometrySettings(max_submissions=1.)
    assert settings.max_submissions == 1 and isinstance(settings.max_submissions, int), (
        'The values should be normalized'
    )
    assert AstrometrySettings().max_submissions is None, 'There should be no limit by default'
    with pytest.raises(PipnickParameterError, match='max_rms'):
        AstrometrySettings(max_rms=0.)
    with pytest.raises(PipnickParameterError, match='wide_search'):
        AstrometrySettings(wide_search=-1.)
