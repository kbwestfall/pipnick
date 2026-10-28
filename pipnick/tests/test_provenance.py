"""
Tests for the provenance hashes.

The hashes are stored in existing products, so they must never change: the
reference values below were computed with the implementations that preceded
:mod:`pipnick.pkg.provenance`.  If one of these tests fails, every product
that records that hash would be considered out of date.
"""

import hashlib

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from pipnick.pipelines import astrometry, calibrate
from pipnick.pkg import provenance


def test_digest():
    assert provenance.digest('abc') == hashlib.sha1(b'abc').hexdigest()[:16], (
        'The digest should be the first 16 hexadecimal digits of the SHA-1 hash'
    )
    assert provenance.digest('ab', b'c') == provenance.digest('abc'), (
        'The digest of several parts should be the digest of their concatenation'
    )


def test_reference_values():
    values = {'exptime': 30.0, 'filter': 'V', 'object': 'M11', 'airmass': None, 'ra': '12:00:00'}
    assert provenance.values_hash(values) == '9728dab81116a8bf', (
        'The hash of the table values (PNTABH) should not change'
    )
    names = ['/a/b/n2001.fits', 'n1002.fits', 'c/n1001.fits']
    assert provenance.names_hash(names) == 'a5f58aa4f3efde5a', (
        'The hash of the raw file names (PNRAWHSH) should not change'
    )
    content = {'b': [1, 2.5, 'x'], 'a': {'z': None, 'y': True}}
    assert provenance.json_hash(content) == '5af1ced572b466df', (
        'The hash of the settings should not change'
    )
    cal = {'calibration': {'fit_hash': 'abc', 'night': '2026-08-19'}, 'bands': {'V': {'zp': 24.3}}}
    assert calibrate.calibration_hash(cal) == 'c5d053bb3715f943', (
        'The hash of a calibration file (fit_hash) should not change'
    )


def test_file_hash(tmp_path):
    path = tmp_path / 'file.txt'
    path.write_bytes(b'pipnick\n')
    assert provenance.file_hash(path) == 'ab12a385b2cf64d4', (
        'The hash of a file should not change'
    )


def test_input_hash():
    data = np.arange(12, dtype=float).reshape(3, 4)
    sources = Table({'x': [1.23456, 5.6789], 'y': [9.87654, 3.21]})
    hints = {'scale_est': 0.333, 'center_ra': 10.5, 'radius': 0.2}
    assert astrometry._input_hash(data, sources, hints, False) == 'c5bbe84e3aa6f836', (
        'The hash of the astrometry.net input (sources) should not change; the cache depends '
        'on it'
    )
    assert astrometry._input_hash(data, sources, hints, True) == 'e60e6ae9ff4bb8b8', (
        'The hash of the astrometry.net input (image) should not change; the cache depends '
        'on it'
    )


def test_wcs_hash():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.crval = [150., 2.]
    wcs.wcs.crpix = [512.5, 512.5]
    wcs.wcs.cd = [[-9.25e-5, 0.], [0., 9.25e-5]]
    expected = hashlib.sha1(wcs.to_header(relax=True).tostring().encode()).hexdigest()[:16]
    assert provenance.wcs_hash(wcs) == expected, (
        'The hash of a WCS should be the digest of its FITS header'
    )


def test_header_hashes():
    hdr = fits.Header(
        {'PNRAWHSH': 'aaaa', 'PNTABH': 'bbbb', 'PNASTMTH': 'gaia', 'PNASTWH': 'cccc'}
    )
    assert provenance.image_hash(hdr) == 'aaaa-bbbb', (
        'The image hash should combine PNRAWHSH and PNTABH'
    )
    assert provenance.image_hash(fits.Header()) == 'none-none', (
        'Missing provenance keywords should be "none"'
    )
    assert provenance.solution_hash(hdr) == 'cccc', 'The solution hash should be PNASTWH'
    assert provenance.solution_hash(fits.Header({'PNASTWH': 'cccc'})) == 'none', (
        'Only a pipnick solution (with PNASTMTH) should count'
    )
