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
    tbl, _ = metadata.build_camera_datatable(tmp_path, path_key='rawdir')

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
    tbl, _ = metadata.build_camera_datatable(tmp_path)
    assert tbl['airmass'].mask[0] and tbl['date'].mask[0] and tbl['mjd'].mask[0]


def test_setups(tmp_path):
    # Two setups, with labels assigned in order of appearance
    write_andor(tmp_path / 'n3001.fits', CDELT1U=1, CDELT2U=1, AMPCONF='0.05MHz')
    write_andor(tmp_path / 'n3002.fits')
    write_andor(tmp_path / 'n3003.fits', CDELT1U=1, CDELT2U=1, AMPCONF='0.05MHz')
    write_andor(tmp_path / 'n3004.fits', WINDOW='960 1205 840 1069')
    tbl, _ = metadata.build_camera_datatable(tmp_path)

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
    tbl, _ = metadata.build_camera_datatable(tmp_path, camera='NickelAndor')
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


def write_andor_image(path, level, **kwargs):
    """
    Write a synthetic 4x4-binned Andor frame with a uniform level.
    """
    hdr = fits.Header({'CAMID': 'headmodel DZ936_BV serial 31382', 'OBSTYPE': '?',
                       'OBJECT': 'Dome flat V', 'XPOSURE': 3.0, 'FILTER': 'V',
                       'DATE-BEG': '2026-08-20T04:12:00.000', 'CDELT1U': 4, 'CDELT2U': 4,
                       'AMPCONF': '1.0MHz', 'WINDOW': '1 2088 1 2048', 'CURRTEMP': -79.9})
    hdr.update(kwargs)
    data = np.full((512, 522), level, dtype=np.uint16)
    fits.PrimaryHDU(data=data, header=hdr).writeto(path)


def test_quality(tmp_path):
    write_andor_image(tmp_path / 'n1.fits', 20000)
    # Too little light for a flat
    write_andor_image(tmp_path / 'n2.fits', 900)
    # Saturated
    write_andor_image(tmp_path / 'n3.fits', 65500)
    # Empty image
    write_andor_image(tmp_path / 'n4.fits', 0)
    # A closed-shutter bias with a light leak
    write_andor_image(tmp_path / 'n5.fits', 5000, OBSTYPE='DARK', OBJECT='bias', XPOSURE=0.2)
    # Labeled as a bias, with the shutter open, but at the bias level
    write_andor_image(tmp_path / 'n6.fits', 301, OBJECT='bias', XPOSURE=0.2)
    # A warm detector
    write_andor_image(tmp_path / 'n7.fits', 20000, CURRTEMP=5.8)

    # The header checks only catch the warm detector
    tbl, rejected = metadata.build_camera_datatable(tmp_path, classify='header')
    assert list(rejected) == [False] * 6 + [True]
    assert 'too warm' in tbl['note'][6]
    assert tbl['frametype'][5] == 'unknown'

    tbl, rejected = metadata.build_camera_datatable(tmp_path, classify='stats')
    assert list(rejected) == [False, True, True, True, True, False, True]
    assert 'too low for a flat' in tbl['note'][1]
    assert 'saturated' in tbl['note'][2]
    assert 'empty' in tbl['note'][3]
    assert 'too high for a bias' in tbl['note'][4]
    assert tbl['frametype'][5] == 'bias' and 'reclassified' in tbl['note'][5]


def test_commented_round_trip(tmp_path):
    write_andor(tmp_path / 'n2001.fits', OBJECT='test bias', FILTER="i'")
    write_andor(tmp_path / 'n2002.fits')
    write_andor(tmp_path / 'n2003.fits', OBJECT='NGC 7331 I', FILTER='I')
    tbl, _ = metadata.build_camera_datatable(tmp_path)
    ofile = tmp_path / 'test.ecsv'
    commented = np.array([True, False, True])
    metadata.write_metadata(tbl, ofile, commented=commented)

    # A comment line that is not a table row is ignored
    with open(ofile, 'a') as f:
        f.write('# n2002 looks a bit odd\n')

    # Commented rows are excluded by default
    assert list(metadata.read_metadata(ofile)['file']) == ['n2002.fits']
    # ... or included, including values with quotes and spaces
    back, back_commented = metadata.read_metadata(ofile, include_commented=True)
    assert list(back['file']) == ['n2001.fits', 'n2002.fits', 'n2003.fits']
    assert np.array_equal(back_commented, commented)
    assert back['filter'][0] == "i'" and back['object'][2] == 'NGC 7331 I'

    # Updating the meta dictionary preserves the commented rows
    metadata.update_metadata_meta(ofile, rdxdir='/tmp/rdx')
    back, back_commented = metadata.read_metadata(ofile, include_commented=True)
    assert back.meta['rdxdir'] == '/tmp/rdx' and np.array_equal(back_commented, commented)


def test_append(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    ofile = tmp_path / 'raw.ecsv'
    write_andor(rawdir / 'n2001.fits')
    write_andor(rawdir / 'n2002.fits', AMPCONF='0.05MHz')
    write_andor(rawdir / 'n2003.fits')
    metadata.build_metadata(rawdir=rawdir, filename=ofile)

    # User edits: comment out one row and change a frame type
    tbl, commented = metadata.read_metadata(ofile, include_commented=True)
    commented[0] = True
    tbl['frametype'][2] = 'test'
    metadata.write_metadata(tbl, ofile, commented=commented)

    # New files, including one with a new setup
    write_andor(rawdir / 'n2004.fits', AMPCONF='0.05MHz')
    write_andor(rawdir / 'n2005.fits', CDELT1U=1, CDELT2U=1)
    _, used = metadata.build_metadata(rawdir=rawdir, filename=ofile, append=True)
    assert list(used['file']) == ['n2002.fits', 'n2003.fits', 'n2004.fits', 'n2005.fits']

    tbl, commented = metadata.read_metadata(ofile, include_commented=True)
    assert list(tbl['file']) == [f'n200{i}.fits' for i in range(1, 6)]
    assert list(commented) == [True, False, False, False, False]
    assert tbl['frametype'][2] == 'test'
    # Existing setups keep their labels
    assert list(tbl['setup']) == ['A', 'B', 'A', 'B', 'C']
    assert tbl.meta['setups']['C']['xbin'] == 1

    # Appending again changes nothing
    metadata.build_metadata(rawdir=rawdir, filename=ofile, append=True)
    again, again_commented = metadata.read_metadata(ofile, include_commented=True)
    assert len(again) == 5 and np.array_equal(again_commented, commented)

    # Overwriting discards the edits
    metadata.build_metadata(rawdir=rawdir, filename=ofile, overwrite=True)
    tbl, commented = metadata.read_metadata(ofile, include_commented=True)
    assert not np.any(commented) and tbl['frametype'][2] == 'science'

    with pytest.raises(ValueError):
        metadata.build_metadata(rawdir=rawdir, filename=ofile, overwrite=True, append=True)
    # Cannot append files from a different directory
    with pytest.raises(ValueError):
        metadata.build_metadata(rawdir=tmp_path, filename=ofile, append=True)


def test_summary(tmp_path):
    write_andor(tmp_path / 'n2001.fits', OBSTYPE='DARK', OBJECT='bias', XPOSURE=0.2)
    write_andor(tmp_path / 'n2002.fits', OBJECT='Dome flat V')
    write_andor(tmp_path / 'n2003.fits', FILTER='V')
    write_andor(tmp_path / 'n2004.fits', FILTER='R')
    write_andor(tmp_path / 'n2005.fits', AMPCONF='0.05MHz')
    tbl, _ = metadata.build_camera_datatable(tmp_path)
    summary = metadata.summarize(tbl, commented=np.array([False, False, False, False, True]))
    assert 'Setup A has R science frames but no R flats.' in summary
    assert 'V flats' not in summary
    # Setup B is commented out entirely
    assert 'Setup B' not in summary.split('Warnings:')[1]
    assert 'Commented rows: 1' in summary
