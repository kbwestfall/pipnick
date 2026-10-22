"""
Tests for the reduction pipeline, using synthetic FITS files.
"""

import numpy as np
import pytest
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.table import Table

from pipnick import cameras, metadata
from pipnick.pipelines import reduction
from pipnick.pipelines.reduction import check_comb_groups, check_single_setup, reduce_all
from pipnick.pkg.exceptions import PipnickDataModelError
from pipnick.tests.test_metadata import write_andor, write_andor_image


def test_one_setup(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    write_andor(rawdir / 'n1001.fits')
    write_andor(rawdir / 'n1002.fits', AMPCONF='0.05MHz')
    table = tmp_path / 'raw.ecsv'
    metadata.build_metadata(rawdir=rawdir, filename=table)

    # Data from more than one setup cannot be reduced together
    with pytest.raises(PipnickDataModelError, match='more than one setup'):
        reduce_all(table, rdxdir=tmp_path / 'rdx')


def test_per_setup_file(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    write_andor(rawdir / 'n1001.fits')
    write_andor(rawdir / 'n1002.fits', AMPCONF='0.05MHz')
    files, _ = metadata.build_metadata(rawdir=rawdir, filename=tmp_path / 'n1.ecsv', split=True)
    assert len(files) == 2, 'Each setup should be written to its own file'

    # Each per-setup file passes the setup check
    for f in files:
        _, tbl = metadata.build_metadata(filename=f)
        check_single_setup(tbl, f)


def test_single_file_table(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    write_andor(rawdir / 'n1001.fits')
    table = tmp_path / 'raw.ecsv'
    metadata.build_metadata(rawdir=rawdir, filename=table)

    # A table with all setups has no combination groups
    with pytest.raises(PipnickDataModelError, match='--split'):
        reduce_all(table, rdxdir=tmp_path / 'rdx')


def comb_table(rows):
    """
    Construct a minimal metadata table for checking the combination groups
    from a list of (frametype, comb_id, filter, exptime) tuples.
    """
    frametype, comb_id, filt, exptime = (list(c) for c in zip(*rows))
    return Table(
        {
            'file': [f'n{i:04d}.fits' for i in range(len(rows))], 'frametype': frametype,
            'comb_id': comb_id, 'filter': filt, 'exptime': exptime,
            'mjd': 61000. + np.arange(len(rows))
        }
    )


valid_groups = [
    ('bias', 1, 'V', 0.2), ('bias', 1, 'V', 0.2), ('domeflat', 2, 'V', 3.),
    ('domeflat', 2, 'V', 3.), ('skyflat', 3, 'V', 10.), ('domeflat', 4, 'R', 3.),
    ('dark', 5, 'V', 120.004), ('dark', 5, 'V', 120.006), ('science', 6, 'V', 30.),
    ('science', 6, 'V', 30.), ('focus', 7, 'V', 5.), ('focus', 7, 'R', 5.)
]


@pytest.mark.parametrize(
    'rows, message', [
        ([('bias', 101, 'V', 0.2), ('dark', 101, 'V', 0.2)], 'mixes frame types'),
        ([('bias', 101, 'V', 0.2)], 'more than one bias group'),
        ([('dark', 101, 'V', 120.), ('dark', 101, 'V', 60.)], 'mixes exposure times'),
        ([('dark', 101, 'V', 120.02)], 'more than one group of 120 s darks'),
        ([('domeflat', 101, 'B', 3.), ('domeflat', 101, 'R', 3.)], 'mixes filters'),
        ([('skyflat', 101, 'V', 5.)], 'more than one group of V skyflat'),
        ([('science', 101, 'V', 30.), ('science', 101, 'R', 30.)], 'mixes filters')
    ]
)
def test_comb_group_violations(rows, message):
    # Each set of rows breaks one rule when added to a valid set of groups
    with pytest.raises(PipnickDataModelError, match=message):
        check_comb_groups(comb_table(valid_groups + rows), 'test.ecsv')


def test_comb_groups(caplog):
    check_comb_groups(comb_table(valid_groups), 'test.ecsv')
    assert 'not reduced' in caplog.text, 'Frames with other types should be noted'

    # Different flat types can share a filter, but only one is used
    check_comb_groups(comb_table([('domeflat', 1, 'V', 3.), ('skyflat', 2, 'V', 3.)]), 'test')

    # Violations are all listed together
    rows = [('bias', 1, 'V', 0.2), ('bias', 2, 'V', 0.2), ('science', 3, 'V', 30.),
            ('science', 3, 'R', 30.)]
    with pytest.raises(PipnickDataModelError) as err:
        check_comb_groups(comb_table(rows), 'test.ecsv')
    assert 'more than one bias group' in str(err.value) and 'mixes filters' in str(err.value), (
        'All violations should be reported'
    )

    # Mixing exposure times in a science group is only a warning
    caplog.clear()
    rows = [('science', 1, 'V', 30.), ('science', 1, 'V', 60.)]
    check_comb_groups(comb_table(rows), 'test.ecsv')
    assert 'mixes exposure times' in caplog.text, (
        'A science group with different exposure times should issue a warning'
    )


def test_names():
    assert reduction._safe_name('NGC 7331 I') == 'NGC_7331_I', 'Spaces should become underscores'
    assert reduction._safe_name("Barnard's") == 'Barnards', 'Quotes should be removed'
    assert reduction._safe_name('SN2026abc+1.5-x/y') == 'SN2026abc+1.5-xy', (
        'Only letters, digits, underscores, periods, plus, and minus should be kept'
    )
    assert reduction._dark_name(120.004) == 'Dark_120s.fits', (
        'Dark exposure times should be rounded to 0.1 s, without trailing zeros'
    )
    assert reduction._dark_name(4.31) == 'Dark_4.3s.fits', (
        'Dark exposure times should be rounded to 0.1 s'
    )
    assert reduction._flat_name('skyflat', 'V') == 'SkyFlat_V.fits', (
        'Flat names should give the type and the filter'
    )
    assert reduction._flat_name('skyflat', "i'") == 'SkyFlat_ip.fits', (
        'Primed filters should be named with a p'
    )
    assert reduction._flat_name('domeflat', 'i') == 'DomeFlat_i.fits', (
        'Unprimed filters should be distinct from primed filters'
    )
    assert reduction._flat_name('domeflat', '') == 'DomeFlat_none.fits', (
        'Flats with an unknown filter should be named with none'
    )

    tbl = Table(
        {
            'object': ['M11 V', 'focus', 'M11 V', 'M 11'], 'filter': ["i'"] * 4,
            'comb_id': [7, 7, 7, 7],
            'date': ['2026-08-21T00:10:00.000', '2026-08-20T23:50:00.000',
                     '2026-08-21T00:20:00.000', '2026-08-20T23:40:00.000'],
            'mjd': [3., 1., 4., 0.]
        }
    )
    indx = reduction._group_indices(tbl)[7]
    assert list(indx) == [3, 1, 0, 2], 'Frames in a group should be ordered by MJD'
    assert reduction._science_name(tbl, indx, 7) == 'M11_V_ip_C007_2026-08-20.fits', (
        'Science names should use the most common object, the filter, and the UT date of the '
        'first frame'
    )
    tbl['object'] = ['b', 'a', 'b', 'a']
    assert reduction._science_name(tbl, indx, 7) == 'a_ip_C007_2026-08-20.fits', (
        'Ties in the most common object should go to the earliest frame'
    )


def write_night(rawdir):
    """
    Write a synthetic night of 4x4-binned Andor frames with biases, V dome
    flats, V sky flats, an R dome flat, a pair of 120 s darks, and two V
    science frames.
    """
    rawdir.mkdir()
    frames = [
        ('n5001', 300, 'bias', 0.2, 'V', 'DARK'), ('n5002', 300, 'bias', 0.2, 'V', 'DARK'),
        ('n5003', 300, 'bias', 0.2, 'V', 'DARK'), ('n5004', 20000, 'Dome flat V', 3., 'V', '?'),
        ('n5005', 20000, 'Dome flat V', 3., 'V', '?'),
        ('n5006', 15000, 'twilight flat', 10., 'V', '?'),
        ('n5007', 15000, 'twilight flat', 10., 'V', '?'),
        ('n5008', 20000, 'Dome flat R', 3., 'R', '?'), ('n5009', 310, 'Dark', 120., 'V', 'DARK'),
        ('n5010', 310, 'Dark', 120., 'V', 'DARK'), ('n5011', 5000, 'M11 V', 30., 'V', '?'),
        ('n5012', 5000, 'M11 V', 30., 'V', '?')
    ]
    for i, (name, level, obj, exptime, filt, obstype) in enumerate(frames):
        write_frame(
            rawdir / f'{name}.fits', level, obj, exptime, filt, obstype, f'04:{i:02d}:00'
        )


def write_frame(path, level, obj, exptime, filt, obstype, time):
    """
    Write a synthetic 4x4-binned Andor frame with coordinates, with a uniform
    level in the data section and a level of 300 ADU in the prescan and
    overscan columns.
    """
    write_andor_image(
        path, level, OBJECT=obj, XPOSURE=exptime, FILTER=filt, OBSTYPE=obstype,
        RA='18:51:05.0', DEC='-06:16:12', **{'DATE-BEG': f'2026-08-20T{time}.000'}
    )
    with fits.open(path, mode='update') as hdu:
        hdu[0].data[:,:5] = 300
        hdu[0].data[:,517:] = 300


def set_comb_id(table, files, comb_id):
    """
    Set the combination group of a set of frames in a per-setup file.
    """
    set_column(table, files, 'comb_id', comb_id)


def set_column(table, files, column, value):
    """
    Set a column value for a set of frames in a per-setup file, as a user
    would by editing the file.
    """
    tbl, commented = metadata.read_metadata(table, include_commented=True, expand=False)
    # Rebuild the column, so that longer strings are not truncated
    tbl[column] = [value if f in files else v for f, v in zip(tbl['file'], tbl[column])]
    metadata.write_metadata(tbl, table, commented=commented)


def comb_id(table, file):
    """
    Return the combination group of a frame in a per-setup file.
    """
    tbl = metadata.read_metadata(table)
    return int(tbl['comb_id'][tbl['file'] == file][0])


def mtimes(rdxdir):
    """
    Return the modification time of each FITS file in a directory.
    """
    return {f.name: f.stat().st_mtime_ns for f in rdxdir.glob('*.fits')}


@pytest.fixture
def night(tmp_path):
    """
    Build the per-setup file for a synthetic night, with the two science frames
    in one group.
    """
    rawdir = tmp_path / 'raw'
    write_night(rawdir)
    files, _ = metadata.build_metadata(rawdir=rawdir, filename=tmp_path / 'n5.ecsv', split=True)
    table = files[0]
    set_comb_id(table, ['n5012.fits'], comb_id(table, 'n5011.fits'))
    return table


def test_reduce(night, tmp_path):
    rdxdir = tmp_path / 'rdx'
    science = reduce_all(night, rdxdir=rdxdir)

    sci_name = f'M11_V_V_C{comb_id(night, "n5011.fits"):03d}_2026-08-20.fits'
    expected = [
        'Bias.fits', 'Dark_120s.fits', 'DomeFlat_R.fits', 'DomeFlat_V.fits', 'SkyFlat_V.fits',
        sci_name
    ]
    assert sorted(f.name for f in rdxdir.glob('*.fits')) == sorted(expected), (
        'One output should be written for each combination group'
    )
    assert [f.name for f in science] == [sci_name], 'The science outputs should be returned'

    bias = fits.getheader(rdxdir / 'Bias.fits')
    assert bias['PNSETUP'] == '4x4_fast_full', 'The setup should be recorded'
    assert bias['PNNRAW'] == 3, 'The number of raw frames should be recorded'
    assert list(bias['HISTORY']) == [
        f'pipnick raw frame: n500{i}.fits' for i in range(1, 4)
    ], 'The raw frames should be listed in the HISTORY'

    flat = fits.getheader(rdxdir / 'SkyFlat_V.fits')
    assert flat['PNBIAS'] == 'Bias.fits' and flat['PNBIASH'] == bias['PNRAWHSH'], (
        'Flats should record the bias they used'
    )

    hdr = fits.getheader(rdxdir / sci_name)
    assert hdr['PNCOMBID'] == comb_id(night, 'n5011.fits'), 'The group should be recorded'
    assert hdr['PNNRAW'] == 2 and hdr['NCOMBINE'] == 2, 'The science frames should be combined'
    assert hdr['PNFLAT'] == 'SkyFlat_V.fits', 'Sky flats should be used over dome flats'
    assert hdr['PNFLATH'] == flat['PNRAWHSH'], 'The hash of the flat should be recorded'
    history = list(hdr['HISTORY'])
    assert (
        'pipnick raw frame: n5011.fits' in history and 'pipnick raw frame: n5012.fits' in history
    ), 'The raw frames should be listed'
    assert 'pipnick calibration: SkyFlat_V.fits' in history, (
        'The calibrations should be listed'
    )
    assert hdr['EXPTIME'] == 30., 'The exposure time should be the median of the group'

    # The bias level is removed and the flat is normalized
    # NOTE: Memory mapping would keep the file open, and Windows does not allow
    # an open file to be overwritten by the reprocessing below
    data = fits.getdata(rdxdir / sci_name, memmap=False)
    assert np.allclose(np.median(data), 4700., rtol=0.01), (
        'The science frames should be bias-subtracted and flat-fielded'
    )

    # Rerunning does not reprocess anything
    before = mtimes(rdxdir)
    reduce_all(night, rdxdir=rdxdir)
    assert mtimes(rdxdir) == before, 'Outputs that are up to date should not be rewritten'

    # Unless overwrite is set
    reduce_all(night, rdxdir=rdxdir, overwrite=True)
    after = mtimes(rdxdir)
    assert all(after[f] != before[f] for f in before), 'Overwrite should reprocess all outputs'


def test_reduce_new_science(night, tmp_path, caplog):
    rdxdir = tmp_path / 'rdx'
    reduce_all(night, rdxdir=rdxdir)
    rawdir = tmp_path / 'raw'

    # A new science frame is added in its own group
    write_frame(rawdir / 'n5013.fits', 5000, 'M11 V', 30., 'V', '?', '05:00:00')
    metadata.build_metadata(rawdir=rawdir, filename=tmp_path / 'n5.ecsv', split=True, append=True)
    before = mtimes(rdxdir)
    science = reduce_all(night, rdxdir=rdxdir)
    new_name = f'M11_V_V_C{comb_id(night, "n5013.fits"):03d}_2026-08-20.fits'
    assert new_name in [f.name for f in science], 'The new frame should be reduced'
    after = mtimes(rdxdir)
    assert all(after[f] == before[f] for f in before), 'Existing outputs should not change'

    # Adding it to the existing science group reprocesses only that group
    sci_name = f'M11_V_V_C{comb_id(night, "n5011.fits"):03d}_2026-08-20.fits'
    set_comb_id(night, ['n5013.fits'], comb_id(night, 'n5011.fits'))
    caplog.clear()
    reduce_all(night, rdxdir=rdxdir)
    final = mtimes(rdxdir)
    changed = sorted(f for f in after if final[f] != after[f])
    assert changed == [sci_name], 'Only the science group with the new frame should be reprocessed'
    assert fits.getheader(rdxdir / sci_name)['PNNRAW'] == 3, 'The new frame should be included'
    assert (rdxdir / new_name).is_file(), 'Stale outputs should not be deleted'
    assert new_name in caplog.text, 'Stale outputs should be reported'


def test_reduce_new_bias(night, tmp_path):
    rdxdir = tmp_path / 'rdx'
    reduce_all(night, rdxdir=rdxdir)
    rawdir = tmp_path / 'raw'

    # A new bias reprocesses everything that depends on it
    write_frame(rawdir / 'n5013.fits', 300, 'bias', 0.2, 'V', 'DARK', '05:00:00')
    metadata.build_metadata(rawdir=rawdir, filename=tmp_path / 'n5.ecsv', split=True, append=True)
    before = mtimes(rdxdir)
    reduce_all(night, rdxdir=rdxdir)
    after = mtimes(rdxdir)
    assert all(after[f] != before[f] for f in before), (
        'A new bias should reprocess the bias and all outputs that use it'
    )
    assert fits.getheader(rdxdir / 'Bias.fits')['PNNRAW'] == 4, 'The new bias should be included'


def test_reduce_other_setup(night, tmp_path):
    rdxdir = tmp_path / 'rdx'
    rdxdir.mkdir()
    hdr = fits.Header({'PNSETUP': '1x1_fast_full'})
    fits.PrimaryHDU(data=np.zeros((4, 4)), header=hdr).writeto(rdxdir / 'Bias.fits')
    with pytest.raises(PipnickDataModelError, match='different setup'):
        reduce_all(night, rdxdir=rdxdir)


def test_frame_stats(tmp_path, caplog):
    path = tmp_path / 'n6001.fits'
    write_frame(path, 5000, 'M11 V', 30., 'V', '?', '04:00:00')
    with fits.open(path, mode='update') as hdu:
        hdu[0].data[100,100:103] = 65535
        hdu[0].data[200,200] = 30000
    ccd = reduction.process_frame(path, flag_cosmics=True)

    assert 'has 3 saturated pixels' in caplog.text, 'Saturated pixels should be counted'
    assert 'mean level is 300.0 ADU' in caplog.text, 'The overscan level should be reported'
    assert ccd.meta['PNOSCAN'] == 300., 'The overscan level should be saved to the header'
    assert f'Trimmed image to {ccd.shape[1]} x {ccd.shape[0]}' in caplog.text, (
        'The size of the trimmed image should be reported'
    )
    assert 'Removed cosmic rays from n6001.fits; 1 pixels affected' in caplog.text, (
        'The number of pixels affected by cosmic rays should be reported'
    )


def test_bias_stats(caplog):
    rng = np.random.default_rng(1)
    gain = 2.
    readnoise = 5.
    frames = [
        CCDData(
            rng.normal(loc=10., scale=readnoise, size=(200, 200)), unit='electron',
            meta={'PROCGAIN': gain, 'PROCRN': readnoise, 'PNOSCAN': 300.}
        ) for _ in range(3)
    ]
    reduction._report_bias_stats(frames, ['b1.fits', 'b2.fits', 'b3.fits'])
    assert 'b1.fits: mean bias level is 305.0 ADU' in caplog.text, (
        'The bias level should include the overscan level and be in ADU'
    )
    measured = float(caplog.text.split('pairs of bias frames: median ')[1].split()[0])
    assert np.isclose(measured, readnoise, rtol=0.05), (
        'The read noise should be the standard deviation of the pair differences divided by '
        'sqrt(2)'
    )


def test_apply_table_values():
    camera = cameras.get_camera('NickelC2')
    ccd = CCDData(
        np.zeros((4, 4)), unit='electron',
        meta=fits.Header({'OBJECT': 'M11 V', 'FILTNAM': 'V', 'EXPTIME': 60, 'RA': '18:51:05.0'})
    )
    values = {'object': 'M11', 'filter': 'R', 'exptime': 60., 'ra': '18:51:05.0', 'dec': None}
    reduction._apply_table_values(ccd, values, camera)

    assert ccd.meta['OBJECT'] == 'M11' and ccd.meta['ORIGOBJ'] == 'M11 V', (
        'A corrected object should be written, keeping the raw value'
    )
    assert ccd.meta['FILTNAM'] == 'R' and ccd.meta['FILTER'] == 'R', (
        'The filter should be written to the camera keyword and the standard keyword'
    )
    assert ccd.meta['ORIGFILT'] == 'V', 'The raw filter should be kept'
    assert 'ORIGEXP' not in ccd.meta and 'ORIGRA' not in ccd.meta, (
        'Raw values that agree with the table should not be duplicated'
    )
    assert 'DEC' not in ccd.meta, 'Unknown table values should not be written'


def test_reduce_table_edits(night, tmp_path, caplog):
    rdxdir = tmp_path / 'rdx'
    reduce_all(night, rdxdir=rdxdir)
    cid = comb_id(night, 'n5011.fits')
    old_name = f'M11_V_V_C{cid:03d}_2026-08-20.fits'
    hdr = fits.getheader(rdxdir / old_name)
    assert hdr['OBJECT'] == 'M11 V' and 'ORIGOBJ' not in hdr, (
        'Unedited table values should match the raw header'
    )
    assert hdr['XPOSURE'] == 30. and hdr['EXPTIME'] == 30., (
        'The exposure time should be written to the camera and standard keywords'
    )

    # Correct the object name in the table
    set_column(night, ['n5011.fits', 'n5012.fits'], 'object', 'M11')
    before = mtimes(rdxdir)
    caplog.clear()
    science = reduce_all(night, rdxdir=rdxdir)
    new_name = f'M11_V_C{cid:03d}_2026-08-20.fits'
    assert [f.name for f in science] == [new_name], 'The corrected object should set the name'
    hdr = fits.getheader(rdxdir / new_name)
    assert hdr['OBJECT'] == 'M11' and hdr['ORIGOBJ'] == 'M11 V', (
        'The corrected object should be in the header, with the raw value kept'
    )
    after = mtimes(rdxdir)
    assert all(after[f] == before[f] for f in before if f != old_name), (
        'Editing the science object should not reprocess the calibrations'
    )
    assert old_name in caplog.text, 'The output with the old name should be reported as stale'

    # Correct the coordinates of one frame; the name does not change
    set_column(night, ['n5011.fits'], 'ra', '18:51:06.0')
    caplog.clear()
    reduce_all(night, rdxdir=rdxdir)
    assert 'metadata table values changed' in caplog.text, (
        'An edited table value should trigger reprocessing'
    )
    hdr = fits.getheader(rdxdir / new_name)
    assert hdr['RA'] == '18:51:06.0' and hdr['ORIGRA'] == '18:51:05.0', (
        'The corrected coordinates should be in the header, with the raw value kept'
    )
    final = mtimes(rdxdir)
    assert [f for f in final if final[f] != after.get(f)] == [new_name], (
        'Only the science output with the edited frame should be reprocessed'
    )
