"""
Tests for the metadata table, using synthetic FITS files.
"""

import numpy as np
import pytest
from astropy.coordinates import EarthLocation
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.utils import iers

from pipnick import cameras
from pipnick.metadata import build, setups
from pipnick.metadata.build import build_camera_datatable, build_metadata, summarize
from pipnick.metadata.groups import mosaic_columns
from pipnick.metadata.setups import find_setup_files, setup_file_keys
from pipnick.metadata.table import (
    metadata_columns, read_metadata, update_metadata_meta, write_metadata
)
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickPathError
from pipnick.pkg.version import version


def write_c2(path, **kwargs):
    """
    Write a small synthetic file from the original Nickel camera.
    """
    hdr = fits.Header(
        {
            'INSTRUME': 'Nickel Direct Camera', 'DSENSOR': 'Loral 2Kx2K', 'OBSTYPE': 'OBJECT',
            'OBJECT': 'NGC_6543', 'EXPTIME': 60, 'FILTNAM': 'R', 'AIRMASS': 1.2,
            'RA': '17:58:33.4', 'DEC': '66:37:59', 'DATE-BEG': '2024-06-27T05:37:01.28',
            'CDELT1U': -2, 'CDELT2U': -2, 'READ-SPD': 80
        }
    )
    hdr.update(kwargs)
    fits.PrimaryHDU(data=np.zeros((4, 4), dtype=np.uint16), header=hdr).writeto(path)


def write_andor(path, **kwargs):
    """
    Write a small synthetic file from the Nickel Andor camera.
    """
    hdr = fits.Header(
        {
            'CAMID': 'headmodel DZ936_BV serial 31382', 'INSTRUME': 'Direct', 'OBSTYPE': '?',
            'OBJECT': 'M11 V', 'XPOSURE': 30.005, 'FILTER': 'V', 'AIRMASS': 1.3,
            'RA': '18:51:05.0', 'DEC': '-06:16:12', 'DATE-BEG': '2026-08-20T04:12:00.000',
            'CDELT1U': 2, 'CDELT2U': 2, 'AMPCONF': '1.0MHz', 'WINDOW': '1 2088 1 2048'
        }
    )
    hdr.update(kwargs)
    fits.PrimaryHDU(data=np.zeros((4, 4), dtype=np.uint16), header=hdr).writeto(path)


def write_andor_image(path, level, **kwargs):
    """
    Write a synthetic 4x4-binned Andor frame with a uniform level.
    """
    hdr = fits.Header(
        {
            'CAMID': 'headmodel DZ936_BV serial 31382', 'OBSTYPE': '?',
            'OBJECT': 'Dome flat V', 'XPOSURE': 3.0, 'FILTER': 'V',
            'DATE-BEG': '2026-08-20T04:12:00.000', 'CDELT1U': 4, 'CDELT2U': 4,
            'AMPCONF': '1.0MHz', 'WINDOW': '1 2088 1 2048', 'CURRTEMP': -79.9
        }
    )
    hdr.update(kwargs)
    data = np.full((512, 522), level, dtype=np.uint16)
    fits.PrimaryHDU(data=data, header=hdr).writeto(path)


def test_columns(tmp_path):
    write_c2(tmp_path / 'd1001.fits', OBJECT='Bias', EXPTIME=0)
    write_c2(tmp_path / 'd1002.fits')
    tbl, _ = build_camera_datatable(tmp_path, path_key='rawdir')

    assert tbl.colnames == list(metadata_columns.keys()), (
        'The table should have the metadata columns, in order'
    )
    assert tbl.meta['camera'] == 'NickelC2', 'The camera should be recorded in the table meta'
    assert tbl.meta['rawdir'] == str(tmp_path), (
        'The data path should be recorded in the table meta'
    )
    assert list(tbl['file']) == ['d1001.fits', 'd1002.fits'], 'Rows should be ordered by file'
    assert list(tbl['frametype']) == ['bias', 'science'], 'The frames should be classified'
    assert tbl['exptime'].dtype.kind == 'f' and tbl['xbin'].dtype.kind == 'i', (
        'The columns should be typed, not strings'
    )
    assert np.all(tbl['xbin'] == 2) and np.all(tbl['ybin'] == 2), (
        'The binning should be positive, regardless of the header sign'
    )
    assert tbl['filter'][1] == 'R', 'The filter should be read from the header'
    assert tbl['mjd'][1] == pytest.approx(60488.23404259), (
        'The MJD should be computed from DATE-BEG'
    )


def test_missing_values(tmp_path):
    write_c2(tmp_path / 'd1001.fits')
    with fits.open(tmp_path / 'd1001.fits', mode='update') as hdu:
        for key in ['AIRMASS', 'DATE-BEG']:
            del hdu[0].header[key]
    tbl, _ = build_camera_datatable(tmp_path)
    assert tbl['airmass'].mask[0] and tbl['date'].mask[0] and tbl['mjd'].mask[0], (
        'Quantities missing from the header should be masked'
    )


def test_computed_airmass(tmp_path):
    write_c2(tmp_path / 'd1001.fits')
    with fits.open(tmp_path / 'd1001.fits', mode='update') as hdu:
        del hdu[0].header['AIRMASS']
    tbl, _ = build_camera_datatable(tmp_path)
    expected = cameras.get_camera('NickelC2').compute_airmass(
        '17:58:33.4', '66:37:59', '2024-06-27T05:37:01.28', exptime=60.
    )
    assert np.isclose(tbl['airmass'][0], expected), (
        'A missing airmass should be computed from the coordinates and the date'
    )
    assert 'airmass computed' in tbl['note'][0], 'A computed airmass should be noted'


def test_compute_airmass():
    camera = cameras.get_camera('NickelC2')
    location = camera.location
    time = Time('2024-06-27T05:00:00', scale='utc')
    with iers.conf.set_temp('auto_download', False):
        lst = time.sidereal_time('apparent', longitude=location.lon)
    zenith = camera.compute_airmass(lst.deg, location.lat.deg, time.isot)
    assert np.isclose(zenith, 1., atol=1e-3), 'The airmass at the zenith should be 1'
    assert camera.compute_airmass('12:00:00', '-80:00:00', time.isot) is None, (
        'The airmass of a target below the horizon should be None'
    )
    assert camera.compute_airmass('not', 'valid', time.isot) is None, (
        'The airmass of invalid coordinates should be None'
    )
    assert camera.compute_airmass('12:00:00', '30:00:00', None) is None, (
        'The airmass without a date should be None'
    )


@pytest.mark.remote_data
def test_observatory_location():
    # The configured location should match astropy's site registry, which
    # requires a download (and so only runs with --remote-data)
    site = EarthLocation.of_site('Lick Observatory')
    for name in cameras.available_cameras():
        location = cameras.get_camera(name).location
        assert location.separation_3d(site).to('m').value < 10., (
            f'The {name} observatory location should match astropy\'s Lick Observatory'
        )


def test_setups(tmp_path):
    # Two setups, with labels assigned in order of appearance
    write_andor(tmp_path / 'n3001.fits', CDELT1U=1, CDELT2U=1, AMPCONF='0.05MHz')
    write_andor(tmp_path / 'n3002.fits')
    write_andor(tmp_path / 'n3003.fits', CDELT1U=1, CDELT2U=1, AMPCONF='0.05MHz')
    write_andor(tmp_path / 'n3004.fits', WINDOW='960 1205 840 1069')
    tbl, _ = build_camera_datatable(tmp_path)

    assert list(tbl['setup']) == ['A', 'B', 'A', 'C'], (
        'Setup labels should be assigned in order of first appearance'
    )
    expected = {
        'A': {
            'camera': 'NickelAndor', 'xbin': 1, 'ybin': 1, 'readmode': 'slow', 'window': 'full',
            'window_name': 'full'
        },
        'B': {
            'camera': 'NickelAndor', 'xbin': 2, 'ybin': 2, 'readmode': 'fast', 'window': 'full',
            'window_name': 'full'
        },
        'C': {
            'camera': 'NickelAndor', 'xbin': 2, 'ybin': 2, 'readmode': 'fast',
            'window': '[960:1205,840:1069]', 'window_name': 'win00'
        }
    }
    assert tbl.meta['setups'] == expected, 'Each setup definition should be recorded'


def test_setup_labels():
    labels = [setups._setup_label(i) for i in [0, 1, 25, 26, 27, 51, 52]]
    assert labels == ['A', 'B', 'Z', 'AA', 'AB', 'AZ', 'BA'], (
        'Setup labels should continue AA, AB, ... after Z'
    )


def test_round_trip(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    write_andor(rawdir / 'n4001.fits', OBSTYPE='DARK', OBJECT='bias', XPOSURE=0.207)
    write_andor(rawdir / 'n4002.fits', AMPCONF='0.05MHz')
    ofile = tmp_path / 'raw.ecsv'
    _, tbl = build_metadata(rawdir=rawdir, filename=ofile)
    back = Table.read(ofile, format='ascii.ecsv')
    assert back.colnames == tbl.colnames, 'The written table should have the same columns'
    assert all(back[c].dtype == tbl[c].dtype for c in tbl.colnames), (
        'The written table should preserve the column types'
    )
    assert back.meta['setups'] == tbl.meta['setups'], (
        'The written table should preserve the setup definitions'
    )
    assert list(back['frametype']) == ['bias', 'science'], (
        'The written table should preserve the frame types'
    )
    # Reloading reads the file instead of rebuilding the table
    _, reloaded = build_metadata(filename=ofile)
    assert list(reloaded['setup']) == ['A', 'B'], 'Reloading should read the existing table'


def test_mixed_cameras(tmp_path):
    write_c2(tmp_path / 'd1001.fits')
    write_andor(tmp_path / 'n1001.fits')
    with pytest.raises(PipnickDataModelError):
        build_camera_datatable(tmp_path)
    tbl, _ = build_camera_datatable(tmp_path, camera='NickelAndor')
    assert list(tbl['file']) == ['n1001.fits'], (
        'Selecting the camera should skip the files from the other camera'
    )


def test_unknown_camera(tmp_path):
    fits.PrimaryHDU(data=np.zeros((4, 4))).writeto(tmp_path / 'x.fits')
    with pytest.raises(PipnickDataModelError):
        build_camera_datatable(tmp_path)


def test_no_files(tmp_path):
    with pytest.raises(PipnickPathError):
        build_camera_datatable(tmp_path)
    with pytest.raises(PipnickPathError):
        build_camera_datatable(tmp_path / 'not_a_directory')
    write_c2(tmp_path / 'd1001.fits')
    with pytest.raises(PipnickPathError):
        build_camera_datatable(tmp_path, camera='NickelAndor')


def test_quality(tmp_path):
    write_andor_image(tmp_path / 'n1.fits', 20000)
    # Too little light for a flat
    write_andor_image(tmp_path / 'n2.fits', 900)
    # Saturated
    write_andor_image(tmp_path / 'n3.fits', 65535)
    # Empty image
    write_andor_image(tmp_path / 'n4.fits', 0)
    # A closed-shutter bias with a light leak
    write_andor_image(tmp_path / 'n5.fits', 5000, OBSTYPE='DARK', OBJECT='bias', XPOSURE=0.2)
    # Labeled as a bias, with the shutter open, but at the bias level
    write_andor_image(tmp_path / 'n6.fits', 301, OBJECT='bias', XPOSURE=0.2)
    # A warm detector
    write_andor_image(tmp_path / 'n7.fits', 20000, CURRTEMP=5.8)

    tbl, rejected = build_camera_datatable(tmp_path, classify='header')
    assert list(rejected) == [False] * 6 + [True], (
        'The header checks should only reject the frame with a warm detector'
    )
    assert 'too warm' in tbl['note'][6], 'The note should explain the rejection'
    assert tbl['frametype'][5] == 'unknown', (
        'The header checks should not reclassify the open-shutter bias'
    )

    tbl, rejected = build_camera_datatable(tmp_path, classify='stats')
    assert list(rejected) == [False, True, True, True, True, False, True], (
        'The image checks should reject the dim flat, the saturated and empty frames, the bias '
        'with a light leak, and the frame with a warm detector'
    )
    assert 'too low for a flat' in tbl['note'][1], 'The note should explain the dim flat'
    assert 'saturated' in tbl['note'][2], 'The note should explain the saturated frame'
    assert 'empty' in tbl['note'][3], 'The note should explain the empty frame'
    assert 'too high for a bias' in tbl['note'][4], 'The note should explain the light leak'
    assert tbl['frametype'][5] == 'bias' and 'reclassified' in tbl['note'][5], (
        'An open-shutter frame labeled as a bias at the bias level should become a bias'
    )


def test_commented_round_trip(tmp_path):
    write_andor(tmp_path / 'n2001.fits', OBJECT='test bias', FILTER="i'")
    write_andor(tmp_path / 'n2002.fits')
    write_andor(tmp_path / 'n2003.fits', OBJECT='NGC 7331 I', FILTER='I')
    tbl, _ = build_camera_datatable(tmp_path)
    ofile = tmp_path / 'test.ecsv'
    commented = np.array([True, False, True])
    write_metadata(tbl, ofile, commented=commented)

    # A comment line that is not a table row is ignored
    with open(ofile, 'a') as f:
        f.write('# n2002 looks a bit odd\n')

    assert list(read_metadata(ofile)['file']) == ['n2002.fits'], (
        'Commented rows should be excluded by default'
    )
    back, back_commented = read_metadata(ofile, include_commented=True)
    assert list(back['file']) == ['n2001.fits', 'n2002.fits', 'n2003.fits'], (
        'Commented rows should be recovered when requested, and other comments ignored'
    )
    assert np.array_equal(back_commented, commented), 'The commented rows should be flagged'
    assert back['filter'][0] == "i'" and back['object'][2] == 'NGC 7331 I', (
        'Values with quotes and spaces should be recovered from commented rows'
    )

    update_metadata_meta(ofile, rdxdir='/tmp/rdx')
    back, back_commented = read_metadata(ofile, include_commented=True)
    assert back.meta['rdxdir'] == '/tmp/rdx', 'The meta dictionary should be updated'
    assert np.array_equal(back_commented, commented), (
        'Updating the meta dictionary should preserve the commented rows'
    )


def test_append(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    ofile = tmp_path / 'raw.ecsv'
    write_andor(rawdir / 'n2001.fits')
    write_andor(rawdir / 'n2002.fits', AMPCONF='0.05MHz')
    write_andor(rawdir / 'n2003.fits')
    build_metadata(rawdir=rawdir, filename=ofile)

    # User edits: comment out one row and change a frame type
    tbl, commented = read_metadata(ofile, include_commented=True)
    commented[0] = True
    tbl['frametype'][2] = 'test'
    write_metadata(tbl, ofile, commented=commented)

    # New files, including one with a new setup
    write_andor(rawdir / 'n2004.fits', AMPCONF='0.05MHz')
    write_andor(rawdir / 'n2005.fits', CDELT1U=1, CDELT2U=1)
    _, used = build_metadata(rawdir=rawdir, filename=ofile, append=True)
    assert list(used['file']) == ['n2002.fits', 'n2003.fits', 'n2004.fits', 'n2005.fits'], (
        'The returned table should exclude the commented rows'
    )

    tbl, commented = read_metadata(ofile, include_commented=True)
    assert list(tbl['file']) == [f'n200{i}.fits' for i in range(1, 6)], (
        'The new files should be appended, ordered by file name'
    )
    assert list(commented) == [True, False, False, False, False], (
        'The commented row should remain commented'
    )
    assert tbl['frametype'][2] == 'test', 'The edited frame type should be preserved'
    assert list(tbl['setup']) == ['A', 'B', 'A', 'B', 'C'], (
        'Existing setups should keep their labels, and the new setup should be labeled next'
    )
    assert tbl.meta['setups']['C']['xbin'] == 1, 'The new setup should be defined'

    build_metadata(rawdir=rawdir, filename=ofile, append=True)
    again, again_commented = read_metadata(ofile, include_commented=True)
    assert len(again) == 5 and np.array_equal(again_commented, commented), (
        'Appending again, with no new files, should change nothing'
    )

    build_metadata(rawdir=rawdir, filename=ofile, overwrite=True)
    tbl, commented = read_metadata(ofile, include_commented=True)
    assert not np.any(commented) and tbl['frametype'][2] == 'science', (
        'Overwriting should rebuild the table, discarding the edits'
    )

    with pytest.raises(ValueError):
        build_metadata(rawdir=rawdir, filename=ofile, overwrite=True, append=True)
    # Cannot append files from a different directory
    with pytest.raises(ValueError):
        build_metadata(rawdir=tmp_path, filename=ofile, append=True)


def test_summary(tmp_path):
    write_andor(tmp_path / 'n2001.fits', OBSTYPE='DARK', OBJECT='bias', XPOSURE=0.2)
    write_andor(tmp_path / 'n2002.fits', OBJECT='Dome flat V')
    write_andor(tmp_path / 'n2003.fits', FILTER='V')
    write_andor(tmp_path / 'n2004.fits', FILTER='R')
    write_andor(tmp_path / 'n2005.fits', AMPCONF='0.05MHz')
    tbl, _ = build_camera_datatable(tmp_path)
    summary = summarize(tbl, commented=np.array([False, False, False, False, True]))
    assert 'Setup A has R science frames but no R flats.' in summary, (
        'Science frames without flats in the same filter should be reported'
    )
    assert 'V flats' not in summary, 'Science frames with flats should not be reported'
    assert 'Setup B' not in summary.split('Warnings:')[1], (
        'Commented-out frames should not be checked for calibrations'
    )
    assert 'Commented rows: 1' in summary, 'The number of commented rows should be reported'


def write_split_night(rawdir):
    """
    Write synthetic Andor files with five setups: two full-frame setups, two
    setups that share a readout window, and one with a second window.
    """
    rawdir.mkdir()
    write_andor(rawdir / 'n3001.fits', CDELT1U=1, CDELT2U=1, AMPCONF='0.05MHz')
    write_andor(rawdir / 'n3002.fits')
    write_andor(rawdir / 'n3003.fits', WINDOW='960 1205 840 1069')
    write_andor(rawdir / 'n3004.fits', WINDOW='960 1205 840 1069', AMPCONF='0.05MHz')
    write_andor(rawdir / 'n3005.fits', WINDOW='101 300 101 300')
    write_andor(rawdir / 'n3006.fits')


def test_split_files(tmp_path):
    rawdir = tmp_path / 'raw'
    write_split_night(rawdir)
    files, used = build_metadata(rawdir=rawdir, filename=tmp_path / 'n3.ecsv', split=True)
    assert all(
        read_metadata(f).meta['pipnick_version'] == version for f in files
    ), 'Each table should record the pipnick version'

    names = sorted(f.name for f in files)
    expected = [
        'n3_1x1_slow_full.ecsv', 'n3_2x2_fast_full.ecsv', 'n3_2x2_fast_win00.ecsv',
        'n3_2x2_fast_win01.ecsv', 'n3_2x2_slow_win00.ecsv'
    ]
    assert names == expected, (
        'Each setup should be written to a file named by its binning, read mode, and window, '
        'and setups with the same window should share its name'
    )
    assert len(used) == 6, 'The returned table should include all setups'
    assert sorted(f.name for f in find_setup_files(tmp_path / 'n3.ecsv')) == expected, (
        'The per-setup files should be found from the base name'
    )

    # The compact layout of a per-setup file
    ofile = tmp_path / 'n3_2x2_fast_full.ecsv'
    compact = read_metadata(ofile, expand=False)
    moved = ['setup'] + setup_file_keys
    expected = [c for c in metadata_columns if c not in moved]
    expected.insert(expected.index('frametype') + 1, 'comb_id')
    for i, name in enumerate(mosaic_columns):
        expected.insert(expected.index('comb_id') + 1 + i, name)
    assert compact.colnames == expected, (
        'The setup values and image size should not be columns of a per-setup file, which '
        'should have the combination groups'
    )
    assert compact.meta['setup'] == {
        'label': 'B', 'camera': 'NickelAndor', 'xbin': 2, 'ybin': 2, 'readmode': 'fast',
        'window': 'full', 'window_name': 'full', 'naxis1': 4, 'naxis2': 4
    }, 'The setup values and image size should be in the file meta dictionary'
    assert 'setups' not in compact.meta, 'A per-setup file should not list all setups'
    assert list(compact['file']) == ['n3002.fits', 'n3006.fits'], (
        'The file should include all frames with its setup'
    )

    # Reading the file restores the columns
    tbl = read_metadata(ofile)
    expected = list(metadata_columns.keys())
    expected.insert(expected.index('setup') + 1, 'comb_id')
    for i, name in enumerate(mosaic_columns):
        expected.insert(expected.index('comb_id') + 1 + i, name)
    assert tbl.colnames == expected, (
        'Reading a per-setup file should restore the columns, in order'
    )
    assert np.all(tbl['mosaic'] == '') and np.all(tbl['wcs_group'] == ''), (
        'By default, no mosaics should be defined'
    )
    assert tbl['xbin'].dtype.kind == 'i' and tbl['readmode'].dtype.kind == 'U', (
        'The restored columns should have their usual types'
    )
    assert list(tbl['setup']) == ['B', 'B'] and list(tbl['naxis1']) == [4, 4], (
        'The restored columns should have the values from the meta dictionary'
    )
    assert list(tbl.meta['setups'].keys()) == ['B'], (
        'The restored meta dictionary should define the single setup'
    )

    # The split files have the same content as a single table
    build_metadata(rawdir=rawdir, filename=tmp_path / 'all.ecsv')
    single = read_metadata(tmp_path / 'all.ecsv')
    combined, _ = build._read_setup_files(files)
    for name in single.colnames:
        assert list(combined[name]) == list(single[name]), (
            f'The {name} column of the per-setup files should match the single table'
        )
    assert combined.meta['setups'] == single.meta['setups'], (
        'The setups in the per-setup files should match the single table'
    )

    # Existing files are not overwritten unless requested
    with pytest.raises(ValueError):
        build_metadata(rawdir=rawdir, filename=tmp_path / 'n3.ecsv', split=True)


def test_split_rejected_setup(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    write_andor(rawdir / 'n1001.fits')
    # The only frame in its setup, rejected because the detector is warm
    write_andor(rawdir / 'n1002.fits', AMPCONF='0.05MHz', CURRTEMP=5.8)
    files, used = build_metadata(rawdir=rawdir, filename=tmp_path / 'n1.ecsv', split=True)
    assert sorted(f.name for f in files) == ['n1_2x2_fast_full.ecsv', 'n1_2x2_slow_full.ecsv'], (
        'A setup with only rejected frames should still be written to a file'
    )
    assert list(used['file']) == ['n1001.fits'], 'The rejected frame should be excluded'
    tbl, commented = read_metadata(
        tmp_path / 'n1_2x2_slow_full.ecsv', include_commented=True
    )
    assert list(tbl['file']) == ['n1002.fits'] and list(commented) == [True], (
        'The rejected frame should be commented out in its file'
    )


def test_split_append(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    base = tmp_path / 'n2.ecsv'
    write_andor(rawdir / 'n2001.fits')
    write_andor(rawdir / 'n2002.fits', AMPCONF='0.05MHz')
    write_andor(rawdir / 'n2003.fits')
    build_metadata(rawdir=rawdir, filename=base, split=True)

    # User edits: comment out one row and change a frame type
    fast = tmp_path / 'n2_2x2_fast_full.ecsv'
    slow = tmp_path / 'n2_2x2_slow_full.ecsv'
    tbl, commented = read_metadata(fast, include_commented=True, expand=False)
    commented[0] = True
    tbl['frametype'][1] = 'test'
    write_metadata(tbl, fast, commented=commented)
    slow_before = slow.read_text()

    # New files, including one with a new setup
    write_andor(rawdir / 'n2004.fits')
    write_andor(rawdir / 'n2005.fits', WINDOW='960 1205 840 1069')
    files, used = build_metadata(rawdir=rawdir, filename=base, split=True, append=True)
    assert sorted(f.name for f in files) == [
        'n2_2x2_fast_full.ecsv', 'n2_2x2_fast_win00.ecsv', 'n2_2x2_slow_full.ecsv'
    ], 'The new setup should be written to a new file'
    assert list(used['file']) == ['n2002.fits', 'n2003.fits', 'n2004.fits', 'n2005.fits'], (
        'The returned table should include all setups, excluding the commented rows'
    )

    tbl, commented = read_metadata(fast, include_commented=True)
    assert list(tbl['file']) == ['n2001.fits', 'n2003.fits', 'n2004.fits'], (
        'The new frame should be added to the file for its setup'
    )
    assert list(commented) == [True, False, False], 'The commented row should stay commented'
    assert tbl['frametype'][1] == 'test', 'The edited frame type should be preserved'
    assert slow.read_text() == slow_before, 'A file without new frames should be unchanged'
    new = read_metadata(tmp_path / 'n2_2x2_fast_win00.ecsv', expand=False)
    assert new.meta['setup']['label'] == 'C', 'The new setup should be labeled after the others'

    # Appending again changes nothing
    contents = {f: f.read_text() for f in files}
    again, _ = build_metadata(rawdir=rawdir, filename=base, split=True, append=True)
    assert sorted(again) == sorted(files), 'Appending again should find the same files'
    assert all(f.read_text() == contents[f] for f in files), (
        'Appending again, with no new files, should change nothing'
    )


def test_split_overwrite(tmp_path, caplog):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    base = tmp_path / 'n2.ecsv'
    write_andor(rawdir / 'n2001.fits')
    write_andor(rawdir / 'n2002.fits', AMPCONF='0.05MHz')
    build_metadata(rawdir=rawdir, filename=base, split=True)

    # Remove the only frame from one setup and rebuild
    (rawdir / 'n2002.fits').unlink()
    files, _ = build_metadata(rawdir=rawdir, filename=base, split=True, overwrite=True)
    assert [f.name for f in files] == ['n2_2x2_fast_full.ecsv'], (
        'Only the remaining setup should be written'
    )
    stale = tmp_path / 'n2_2x2_slow_full.ecsv'
    assert stale.is_file(), 'Existing per-setup files should not be deleted'
    assert stale.name in caplog.text, 'Per-setup files that were not rewritten should be reported'


def test_split_update_meta(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    write_andor(rawdir / 'n2001.fits')
    files, _ = build_metadata(rawdir=rawdir, filename=tmp_path / 'n2.ecsv', split=True)
    update_metadata_meta(files[0], rdxdir='/tmp/rdx')
    compact = read_metadata(files[0], expand=False)
    assert compact.meta['rdxdir'] == '/tmp/rdx', 'The meta dictionary should be updated'
    assert 'xbin' not in compact.colnames and 'setup' in compact.meta, (
        'Updating the meta dictionary should preserve the per-setup layout'
    )


def write_grouping_night(rawdir):
    """
    Write synthetic Andor files with a mix of frame types in one setup.
    """
    rawdir.mkdir()
    for i in range(1, 4):
        write_andor(rawdir / f'n400{i}.fits', OBSTYPE='DARK', OBJECT='bias', XPOSURE=0.2)
    write_andor(rawdir / 'n4004.fits', OBJECT='Dome flat V', FILTER='V', XPOSURE=3.0)
    write_andor(rawdir / 'n4005.fits', OBJECT='Dome flat V', FILTER='V', XPOSURE=3.0)
    write_andor(rawdir / 'n4006.fits', OBJECT='Dome flat R', FILTER='R', XPOSURE=3.0)
    write_andor(rawdir / 'n4007.fits', OBJECT='twilight flat', FILTER='V', XPOSURE=10.0)
    write_andor(rawdir / 'n4008.fits', OBSTYPE='DARK', OBJECT='Dark', XPOSURE=120.004)
    write_andor(rawdir / 'n4009.fits', OBSTYPE='DARK', OBJECT='Dark', XPOSURE=120.006)
    write_andor(rawdir / 'n4010.fits', OBSTYPE='DARK', OBJECT='Dark', XPOSURE=30.005)
    write_andor(rawdir / 'n4011.fits', OBJECT='M11 V')
    write_andor(rawdir / 'n4012.fits', OBJECT='M11 V')
    write_andor(rawdir / 'n4013.fits', OBJECT='Pointing 18')
    write_andor(rawdir / 'n4014.fits', OBJECT='Focus 18')


def test_comb_ids(tmp_path):
    rawdir = tmp_path / 'raw'
    write_grouping_night(rawdir)
    files, used = build_metadata(rawdir=rawdir, filename=tmp_path / 'n4.ecsv', split=True)
    tbl, commented = read_metadata(files[0], include_commented=True)
    ids = dict(zip(tbl['file'], tbl['comb_id']))

    assert tbl['comb_id'].dtype.kind == 'i', 'The combination groups should be integers'
    assert ids['n4001.fits'] == ids['n4002.fits'] == ids['n4003.fits'], (
        'All biases should be combined'
    )
    assert ids['n4004.fits'] == ids['n4005.fits'] != ids['n4006.fits'], (
        'Dome flats should be combined by filter'
    )
    assert ids['n4007.fits'] not in [ids['n4004.fits'], ids['n4006.fits']], (
        'Sky flats should not be combined with dome flats'
    )
    assert ids['n4008.fits'] == ids['n4009.fits'] != ids['n4010.fits'], (
        'Darks should be combined if their exposure times agree to within 0.1 s'
    )
    singles = ['n4011.fits', 'n4012.fits', 'n4013.fits', 'n4014.fits']
    assert len(set(ids[f] for f in singles)) == 4, (
        'Science, pointing, and focus frames should be processed independently'
    )
    assert len(set(tbl['comb_id'])) == 10, (
        'There should be 10 combination groups: biases, V and R dome flats, the sky flat, two '
        'dark exposure times, two science frames, and the pointing and focus frames'
    )

    excluded = list(tbl['file'][commented])
    assert excluded == ['n4013.fits', 'n4014.fits'], (
        'Pointing and focus frames should be commented out by default'
    )
    assert 'excluded by default' in tbl['note'][tbl['file'] == 'n4013.fits'][0], (
        'The note should explain why the pointing frame is commented out'
    )
    assert 'n4013.fits' not in list(used['file']), (
        'Pointing frames should be excluded from the returned table'
    )

    # Neither applies when writing a single table
    _, single = build_metadata(rawdir=rawdir, filename=tmp_path / 'all.ecsv')
    assert 'comb_id' not in single.colnames, 'A single table should not have combination groups'
    assert 'n4013.fits' in list(single['file']), (
        'Pointing frames should not be excluded from a single table'
    )


def test_comb_ids_append(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    base = tmp_path / 'n4.ecsv'
    write_andor(rawdir / 'n4001.fits', OBSTYPE='DARK', OBJECT='bias', XPOSURE=0.2)
    write_andor(rawdir / 'n4002.fits', OBSTYPE='DARK', OBJECT='bias', XPOSURE=0.2)
    write_andor(rawdir / 'n4003.fits', OBJECT='Dome flat V', FILTER='V', XPOSURE=3.0)
    write_andor(rawdir / 'n4004.fits', OBJECT='Dome flat V', FILTER='V', XPOSURE=3.0)
    write_andor(rawdir / 'n4005.fits', OBJECT='Pointing 18')
    files, _ = build_metadata(rawdir=rawdir, filename=base, split=True)

    # User edits: split the V flats into two groups, and restore the pointing frame
    tbl, commented = read_metadata(files[0], include_commented=True, expand=False)
    tbl['comb_id'][tbl['file'] == 'n4004.fits'] = 99
    commented[tbl['file'] == 'n4005.fits'] = False
    write_metadata(tbl, files[0], commented=commented)

    # New frames
    write_andor(rawdir / 'n4006.fits', OBSTYPE='DARK', OBJECT='bias', XPOSURE=0.2)
    write_andor(rawdir / 'n4007.fits', OBJECT='Dome flat V', FILTER='V', XPOSURE=3.0)
    write_andor(rawdir / 'n4008.fits', OBJECT='Dome flat R', FILTER='R', XPOSURE=3.0)
    write_andor(rawdir / 'n4009.fits', OBJECT='Focus 18')
    build_metadata(rawdir=rawdir, filename=base, split=True, append=True)

    tbl, commented = read_metadata(files[0], include_commented=True)
    ids = dict(zip(tbl['file'], tbl['comb_id']))
    flags = dict(zip(tbl['file'], commented))
    assert ids['n4006.fits'] == ids['n4001.fits'], 'A new bias should join the bias group'
    assert ids['n4004.fits'] == 99, 'The edited combination group should be preserved'
    assert ids['n4007.fits'] not in [ids['n4003.fits'], 99], (
        'A new flat should not join a default group that the user split'
    )
    assert ids['n4008.fits'] not in [ids[f] for f in ['n4001.fits', 'n4003.fits', 'n4007.fits']], (
        'A flat with a new filter should get a new group'
    )
    assert not np.ma.is_masked(tbl['comb_id']), 'Every frame should have a group'
    assert not flags['n4005.fits'], 'The restored pointing frame should stay restored'
    assert flags['n4009.fits'], 'A new focus frame should be commented out by default'
