"""
Tests for the names of the pipnick products.
"""

from pathlib import Path

import pytest

from pipnick.pkg import paths
from pipnick.pkg.exceptions import PipnickParameterError


def test_step_dir():
    assert paths.step_dir('rdx', 'photometry') == Path('rdx/photometry'), (
        'The products of a step should be in a subdirectory named for the step'
    )
    with pytest.raises(PipnickParameterError):
        paths.step_dir('rdx', 'unknown')


def test_image_products():
    rdxdir = Path('/data/rdx')
    image = rdxdir / 'M11_V_C014_2026-08-20.fits'
    expected = rdxdir / 'photometry' / 'M11_V_C014_2026-08-20_phot.ecsv'
    assert paths.photometry_table(image) == expected, (
        'The source table should be in photometry/, named for the image'
    )
    expected = rdxdir / 'photometry' / 'M11_V_C014_2026-08-20_removed.ecsv'
    assert paths.removed_table(image) == expected, (
        'The removed sources should be in photometry/, named for the image'
    )
    expected = rdxdir / 'astrometry' / 'M11_V_C014_2026-08-20_points.fits'
    assert paths.points_file(str(image)) == expected, (
        'The user solution should be in astrometry/, named for the image'
    )
    expected = rdxdir / 'calibrate' / 'M11_V_C014_2026-08-20_cal.ecsv'
    assert paths.calibrated_table(image) == expected, (
        'The calibrated table should be in calibrate/, named for the image'
    )


def test_calibration_files():
    expected = Path('rdx/calibrate/2x2_slow_full_2026-08-20.toml')
    assert paths.calibration_file('rdx', '2x2_slow_full', '2026-08-20') == expected, (
        'The calibration file should be named for the setup and night'
    )
    expected = Path('rdx/calibrate/2x2_slow_full_2026-08-20_residuals.ecsv')
    assert paths.residuals_file('rdx', '2x2_slow_full', '2026-08-20') == expected, (
        'The residuals table should be named for the setup and night'
    )


def test_qa_figure():
    expected = Path('rdx/photometry/qa/M11_V_C014_2026-08-20_psf.png')
    assert paths.qa_figure('rdx/photometry', 'M11_V_C014_2026-08-20_psf') == expected, (
        'A QA figure should be in the qa subdirectory of its step, named for its product'
    )
