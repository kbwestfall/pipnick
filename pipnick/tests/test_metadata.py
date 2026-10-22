"""
Tests for the metadata table, using synthetic FITS files.
"""

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from pipnick import metadata
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickPathError


def write_c2(path, **kwargs):
    """
    Write a small synthetic file from the original Nickel camera.
    """
    hdr = fits.Header({'INSTRUME': 'Nickel Direct Camera', 'DSENSOR': 'Loral 2Kx2K',
                       'OBSTYPE': 'OBJECT', 'OBJECT': 'NGC_6543', 'EXPTIME': 60,
                       'FILTNAM': 'R', 'AIRMASS': 1.2, 'RA': '17:58:33.4', 'DEC': '66:37:59',
                       'DATE-BEG': '2024-06-27T05:37:01.28', 'CDELT1U': -2, 'CDELT2U': -2,
                       'READ-SPD': 80})
    hdr.update(kwargs)
    fits.PrimaryHDU(data=np.zeros((4, 4), dtype=np.uint16), header=hdr).writeto(path)


def write_andor(path, **kwargs):
    """
    Write a small synthetic file from the Nickel Andor camera.
    """
    hdr = fits.Header({'CAMID': 'headmodel DZ936_BV serial 31382', 'INSTRUME': 'Direct',
                       'OBSTYPE': '?', 'OBJECT': 'M11 V', 'XPOSURE': 30.005, 'FILTER': 'V',
                       'AIRMASS': 1.3, 'RA': '18:51:05.0', 'DEC': '-06:16:12',
                       'DATE-BEG': '2026-08-20T04:12:00.000', 'CDELT1U': 2, 'CDELT2U': 2,
                       'AMPCONF': '1.0MHz', 'WINDOW': '1 2088 1 2048'})
    hdr.update(kwargs)
    fits.PrimaryHDU(data=np.zeros((4, 4), dtype=np.uint16), header=hdr).writeto(path)


def test_columns(tmp_path):
    write_c2(tmp_path / 'd1001.fits', OBJECT='Bias', EXPTIME=0)
    write_c2(tmp_path / 'd1002.fits')
    tbl = metadata.build_camera_datatable(tmp_path, path_key='rawdir')

    assert tbl.colnames == list(metadata.metadata_columns.keys())
    assert tbl.meta['camera'] == 'NickelC2'
    assert tbl.meta['rawdir'] == str(tmp_path)
    assert list(tbl['file']) == ['d1001.fits', 'd1002.fits']
    assert list(tbl['frametype']) == ['bias', 'science']
    assert tbl['exptime'].dtype.kind == 'f' and tbl['xbin'].dtype.kind == 'i'
    # Binning is positive, regardless of the header sign
    assert np.all(tbl['xbin'] == 2) and np.all(tbl['ybin'] == 2)
    assert tbl['filter'][1] == 'R'
    assert tbl['mjd'][1] == pytest.approx(60488.23404259)


def test_missing_values(tmp_path):
    hdr = {'AIRMASS': None, 'DATE-BEG': None}
    write_c2(tmp_path / 'd1001.fits')
    with fits.open(tmp_path / 'd1001.fits', mode='update') as hdu:
        for key in hdr:
            del hdu[0].header[key]
    tbl = metadata.build_camera_datatable(tmp_path)
    assert tbl['airmass'].mask[0] and tbl['date'].mask[0] and tbl['mjd'].mask[0]


def test_setups(tmp_path):
    # Two setups, with labels assigned in order of appearance
    write_andor(tmp_path / 'n3001.fits', CDELT1U=1, CDELT2U=1, AMPCONF='0.05MHz')
    write_andor(tmp_path / 'n3002.fits')
    write_andor(tmp_path / 'n3003.fits', CDELT1U=1, CDELT2U=1, AMPCONF='0.05MHz')
    write_andor(tmp_path / 'n3004.fits', WINDOW='960 1205 840 1069')
    tbl = metadata.build_camera_datatable(tmp_path)

    assert list(tbl['setup']) == ['A', 'B', 'A', 'C']
    assert tbl.meta['setups'] == {
        'A': {'camera': 'NickelAndor', 'xbin': 1, 'ybin': 1, 'readmode': 'slow',
              'window': 'full'},
        'B': {'camera': 'NickelAndor', 'xbin': 2, 'ybin': 2, 'readmode': 'fast',
              'window': 'full'},
        'C': {'camera': 'NickelAndor', 'xbin': 2, 'ybin': 2, 'readmode': 'fast',
              'window': '[960:1205,840:1069]'},
    }


def test_setup_labels():
    assert [metadata._setup_label(i) for i in [0, 1, 25, 26, 27, 51, 52]] \
                == ['A', 'B', 'Z', 'AA', 'AB', 'AZ', 'BA']


def test_round_trip(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    write_andor(rawdir / 'n4001.fits', OBSTYPE='DARK', OBJECT='bias', XPOSURE=0.207)
    write_andor(rawdir / 'n4002.fits', AMPCONF='0.05MHz')
    ofile = tmp_path / 'raw.ecsv'
    _, tbl = metadata.build_metadata(rawdir=rawdir, filename=ofile)
    back = Table.read(ofile, format='ascii.ecsv')
    assert back.colnames == tbl.colnames
    assert all(back[c].dtype == tbl[c].dtype for c in tbl.colnames)
    assert back.meta['setups'] == tbl.meta['setups']
    assert list(back['frametype']) == ['bias', 'science']
    # Reloading reads the file instead of rebuilding the table
    _, reloaded = metadata.build_metadata(filename=ofile)
    assert list(reloaded['setup']) == ['A', 'B']


def test_mixed_cameras(tmp_path):
    write_c2(tmp_path / 'd1001.fits')
    write_andor(tmp_path / 'n1001.fits')
    with pytest.raises(PipnickDataModelError):
        metadata.build_camera_datatable(tmp_path)
    # Selecting the camera skips the files from the other camera
    tbl = metadata.build_camera_datatable(tmp_path, camera='NickelAndor')
    assert list(tbl['file']) == ['n1001.fits']


def test_unknown_camera(tmp_path):
    fits.PrimaryHDU(data=np.zeros((4, 4))).writeto(tmp_path / 'x.fits')
    with pytest.raises(PipnickDataModelError):
        metadata.build_camera_datatable(tmp_path)


def test_no_files(tmp_path):
    with pytest.raises(PipnickPathError):
        metadata.build_camera_datatable(tmp_path)
    with pytest.raises(PipnickPathError):
        metadata.build_camera_datatable(tmp_path / 'not_a_directory')
    write_c2(tmp_path / 'd1001.fits')
    with pytest.raises(PipnickPathError):
        metadata.build_camera_datatable(tmp_path, camera='NickelAndor')
