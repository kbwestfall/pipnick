"""
Tests for the reduction pipeline, using synthetic FITS files.
"""

import matplotlib.image
import numpy as np
import pytest
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.table import Table

from pipnick import cameras
from pipnick.metadata.build import build_metadata
from pipnick.metadata.names import safe_name
from pipnick.metadata.table import read_metadata, write_metadata
from pipnick.datamodel import dataset
from pipnick.datamodel.dataset import DataSet, check_comb_groups, check_single_setup
from pipnick.pipelines import reduction
from pipnick.reduction import processing
from pipnick.pipelines.reduction import reduce_all
from pipnick.pkg.exceptions import PipnickDataModelError
from pipnick.pkg.report import StepReport
from pipnick.pkg.version import version
from pipnick.tests.test_metadata import write_andor, write_andor_image


def test_one_setup(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    write_andor(rawdir / 'n1001.fits')
    write_andor(rawdir / 'n1002.fits', AMPCONF='0.05MHz')
    table = tmp_path / 'raw.ecsv'
    build_metadata(rawdir=rawdir, filename=table)

    # Data from more than one setup cannot be reduced together
    with pytest.raises(PipnickDataModelError, match='more than one setup'):
        reduce_all(table, rdxdir=tmp_path / 'rdx')


def test_per_setup_file(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    write_andor(rawdir / 'n1001.fits')
    write_andor(rawdir / 'n1002.fits', AMPCONF='0.05MHz')
    files, _ = build_metadata(rawdir=rawdir, filename=tmp_path / 'n1.ecsv', split=True)
    assert len(files) == 2, 'Each setup should be written to its own file'

    # Each per-setup file passes the setup check
    for f in files:
        _, tbl = build_metadata(filename=f)
        check_single_setup(tbl, f)


def test_single_file_table(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    write_andor(rawdir / 'n1001.fits')
    table = tmp_path / 'raw.ecsv'
    build_metadata(rawdir=rawdir, filename=table)

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
    rows = [
        ('bias', 1, 'V', 0.2), ('bias', 2, 'V', 0.2), ('science', 3, 'V', 30.),
        ('science', 3, 'R', 30.)
    ]
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
    assert safe_name('NGC 7331 I') == 'NGC_7331_I', 'Spaces should become underscores'
    assert safe_name("Barnard's") == 'Barnards', 'Quotes should be removed'
    assert safe_name('SN2026abc+1.5-x/y') == 'SN2026abc+1.5-xy', (
        'Only letters, digits, underscores, periods, plus, and minus should be kept'
    )
    assert dataset._dark_name(120.004) == 'Dark_120s.fits', (
        'Dark exposure times should be rounded to 0.1 s, without trailing zeros'
    )
    assert dataset._dark_name(4.31) == 'Dark_4.3s.fits', (
        'Dark exposure times should be rounded to 0.1 s'
    )
    assert dataset._flat_name('skyflat', 'V') == 'SkyFlat_V.fits', (
        'Flat names should give the type and the filter'
    )
    assert dataset._flat_name('skyflat', "i'") == 'SkyFlat_ip.fits', (
        'Primed filters should be named with a p'
    )
    assert dataset._flat_name('domeflat', 'i') == 'DomeFlat_i.fits', (
        'Unprimed filters should be distinct from primed filters'
    )
    assert dataset._flat_name('domeflat', '') == 'DomeFlat_none.fits', (
        'Flats with an unknown filter should be named with none'
    )

    tbl = Table(
        {
            'object': ['M11 V', 'focus', 'M11 V', 'M 11'], 'filter': ["i'"] * 4,
            'comb_id': [7, 7, 7, 7],
            'date': [
                '2026-08-21T00:10:00.000', '2026-08-20T23:50:00.000', '2026-08-21T00:20:00.000',
                '2026-08-20T23:40:00.000'
            ],
            'mjd': [3., 1., 4., 0.]
        }
    )
    indx = dataset._group_indices(tbl)[7]
    assert list(indx) == [3, 1, 0, 2], 'Frames in a group should be ordered by MJD'
    assert dataset._science_name(tbl, indx, 7) == 'M11_V_ip_C007_2026-08-20.fits', (
        'Science names should use the most common object, the filter, and the UT date of the '
        'first frame'
    )
    tbl['object'] = ['b', 'a', 'b', 'a']
    assert dataset._science_name(tbl, indx, 7) == 'a_ip_C007_2026-08-20.fits', (
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
    tbl, commented = read_metadata(table, include_commented=True, expand=False)
    # Rebuild the column, so that longer strings are not truncated
    tbl[column] = [value if f in files else v for f, v in zip(tbl['file'], tbl[column])]
    write_metadata(tbl, table, commented=commented)


def comb_id(table, file):
    """
    Return the combination group of a frame in a per-setup file.
    """
    tbl = read_metadata(table)
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
    files, _ = build_metadata(rawdir=rawdir, filename=tmp_path / 'n5.ecsv', split=True)
    table = files[0]
    set_comb_id(table, ['n5012.fits'], comb_id(table, 'n5011.fits'))
    return table


def test_reduce(night, tmp_path):
    rdxdir = tmp_path / 'rdx'
    science = reduce_all(night, rdxdir=rdxdir, qa=False)

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
    gain = cameras.get_camera('NickelAndor').readout(
        fits.getheader(tmp_path / 'raw' / 'n5011.fits')
    ).gain
    assert np.allclose(np.median(data), 4700. * gain, rtol=0.01), (
        'The science frames should be bias-subtracted, flat-fielded, and in electrons'
    )

    # Rerunning does not reprocess anything
    before = mtimes(rdxdir)
    reduce_all(night, rdxdir=rdxdir, qa=False)
    assert mtimes(rdxdir) == before, 'Outputs that are up to date should not be rewritten'

    # Unless overwrite is set
    reduce_all(night, rdxdir=rdxdir, overwrite=True, qa=False)
    after = mtimes(rdxdir)
    assert all(after[f] != before[f] for f in before), 'Overwrite should reprocess all outputs'


def test_reduce_new_science(night, tmp_path, caplog):
    rdxdir = tmp_path / 'rdx'
    reduce_all(night, rdxdir=rdxdir, qa=False)
    rawdir = tmp_path / 'raw'

    # A new science frame is added in its own group
    write_frame(rawdir / 'n5013.fits', 5000, 'M11 V', 30., 'V', '?', '05:00:00')
    build_metadata(rawdir=rawdir, filename=tmp_path / 'n5.ecsv', split=True, append=True)
    before = mtimes(rdxdir)
    science = reduce_all(night, rdxdir=rdxdir, qa=False)
    new_name = f'M11_V_V_C{comb_id(night, "n5013.fits"):03d}_2026-08-20.fits'
    assert new_name in [f.name for f in science], 'The new frame should be reduced'
    after = mtimes(rdxdir)
    assert all(after[f] == before[f] for f in before), 'Existing outputs should not change'

    # Adding it to the existing science group reprocesses only that group
    sci_name = f'M11_V_V_C{comb_id(night, "n5011.fits"):03d}_2026-08-20.fits'
    set_comb_id(night, ['n5013.fits'], comb_id(night, 'n5011.fits'))
    caplog.clear()
    reduce_all(night, rdxdir=rdxdir, qa=False)
    final = mtimes(rdxdir)
    changed = sorted(f for f in after if final[f] != after[f])
    assert changed == [sci_name], 'Only the science group with the new frame should be reprocessed'
    assert fits.getheader(rdxdir / sci_name)['PNNRAW'] == 3, 'The new frame should be included'
    assert (rdxdir / new_name).is_file(), 'Stale outputs should not be deleted'
    assert new_name in caplog.text, 'Stale outputs should be reported'


def test_reduce_new_bias(night, tmp_path):
    rdxdir = tmp_path / 'rdx'
    reduce_all(night, rdxdir=rdxdir, qa=False)
    rawdir = tmp_path / 'raw'

    # A new bias reprocesses everything that depends on it
    write_frame(rawdir / 'n5013.fits', 300, 'bias', 0.2, 'V', 'DARK', '05:00:00')
    build_metadata(rawdir=rawdir, filename=tmp_path / 'n5.ecsv', split=True, append=True)
    before = mtimes(rdxdir)
    reduce_all(night, rdxdir=rdxdir, qa=False)
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
    ccd = processing.process_frame(path, flag_cosmics=True)

    assert ccd.meta['PNNSAT'] == 3, 'Saturated pixels should be counted'
    assert ccd.meta['PNOSCAN'] == 300., 'The overscan level should be saved to the header'
    assert ccd.meta['PNNCR'] == 1, 'The pixels affected by cosmic rays should be counted'
    assert 'n6001.fits: overscan 300.0 ADU; 3 saturated pixels; 1 cosmic-ray pixels' in (
        caplog.text
    ), 'One line should summarize the processing of the frame, naming it'


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
    level, measured = processing.report_bias_stats(frames, ['b1.fits', 'b2.fits', 'b3.fits'])
    assert np.isclose(level, 305., atol=0.1), (
        'The bias level should include the overscan level and be in ADU'
    )
    assert 'Bias level of the 3 frames: median 305.0 ADU' in caplog.text, (
        'The bias level should be reported'
    )
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
    reduce_all(night, rdxdir=rdxdir, qa=False)
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
    science = reduce_all(night, rdxdir=rdxdir, qa=False)
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
    reduce_all(night, rdxdir=rdxdir, qa=False)
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


def test_reduce_airmass(night, tmp_path):
    rdxdir = tmp_path / 'rdx'
    reduce_all(night, rdxdir=rdxdir, qa=False)
    cid = comb_id(night, 'n5011.fits')
    name = f'M11_V_V_C{cid:03d}_2026-08-20.fits'
    tbl = read_metadata(night)
    airmass = tbl['airmass'][np.isin(tbl['file'], ['n5011.fits', 'n5012.fits'])]
    assert np.isclose(fits.getheader(rdxdir / name)['AIRMASS'], np.mean(airmass)), (
        'The airmass of a combined image should be the mean of its frames'
    )

    # An edited airmass is propagated
    set_column(night, ['n5011.fits', 'n5012.fits'], 'airmass', 1.5)
    before = mtimes(rdxdir)
    reduce_all(night, rdxdir=rdxdir, qa=False)
    after = mtimes(rdxdir)
    assert [f for f in after if after[f] != before[f]] == [name], (
        'Editing the airmass should reprocess only the science image'
    )
    hdr = fits.getheader(rdxdir / name)
    assert hdr['AIRMASS'] == 1.5, 'The edited airmass should be in the header'
    # The synthetic raw frames have no airmass (it is computed for the
    # metadata table), so there is no original value to keep
    assert 'ORIGAIR' not in hdr, 'There is no raw header airmass to keep'
    assert 'airmass computed' in tbl['note'][tbl['file'] == 'n5011.fits'][0], (
        'The airmass of the synthetic frames should be computed'
    )


def test_uncertainty(night, tmp_path):
    rdxdir = tmp_path / 'rdx'
    reduce_all(night, rdxdir=rdxdir, qa=False)
    cid = comb_id(night, 'n5011.fits')

    # The Andor gain and read noise, and the number of binned overscan columns
    # averaged in each row
    raw = tmp_path / 'raw' / 'n5011.fits'
    camera = cameras.get_camera('NickelAndor')
    readout = camera.readout(fits.getheader(raw))
    readnoise = readout.readnoise
    _, biassec = camera.sections(fits.getheader(raw))
    x1, x2 = (int(v) for v in biassec.strip('[]').split(',')[0].split(':'))
    var_rn = readnoise**2 * (1 + 1 / (x2 - x1 + 1))

    # The synthetic frames have no noise, so the uncertainties follow exactly
    # from the read noise and the counts (in electrons): 3 biases at the
    # overscan level; 2 sky flats with 14700 ADU; and 2 science frames with
    # 4700 ADU
    flat = 14700. * readout.gain
    science = 4700. * readout.gain
    var_bias = var_rn / 3
    var_flat = (var_rn + var_bias + flat) / 2
    var_frame = var_rn + var_bias + science + science**2 * var_flat / flat**2
    expected = {
        'Bias.fits': np.sqrt(var_bias),
        f'M11_V_V_C{cid:03d}_2026-08-20.fits': np.sqrt(var_frame / 2)
    }
    for name, sigma in expected.items():
        ccd = CCDData.read(rdxdir / name)
        assert np.allclose(np.median(ccd.uncertainty.array), sigma, rtol=1e-3), (
            f'The uncertainty of {name} should be propagated from the read noise, the overscan, '
            'the counts, and the calibrations'
        )


def test_reduce_summary_and_verdict(night, tmp_path, caplog):
    rdxdir = tmp_path / 'rdx'
    report = StepReport('REDUCTION')
    reduce_all(night, rdxdir=rdxdir, report=report, qa=False)

    summary = Table.read(rdxdir / 'reduce' / 'summary.ecsv')
    assert list(summary['status']) == ['built'] * 6, 'Every product should be built'
    assert summary.meta['pipnick_version'] == version, 'The summary should record the version'
    bias = summary[summary['frametype'] == 'bias'][0]
    assert bias['bias_level'] == pytest.approx(300., abs=0.5), (
        'The bias level should be measured in ADU'
    )
    assert bias['read_noise'] < 1., 'The synthetic biases have no noise'
    science = summary[summary['frametype'] == 'science'][0]
    assert science['flat'] == 'SkyFlat_V.fits' and science['bias'] == 'Bias.fits', (
        'The summary should list the calibrations used'
    )
    assert science['sky'] > 0 and science['cr_frac'] >= 0, (
        'The sky level and cosmic-ray fraction should be measured'
    )
    hdr = fits.getheader(rdxdir / 'Bias.fits')
    assert hdr['PNVERS'] == version, 'The products should record the pipnick version'
    assert hdr['PNBIASLV'] == pytest.approx(bias['bias_level']), (
        'The bias level should be saved to the header of the bias'
    )

    # The noiseless biases have a read noise far below the nominal value
    assert report.status == 'NEEDS ATTENTION' and report.exit_status == 0, (
        'A measurement outside the expected range should need attention, not fail'
    )
    assert any('read noise' in p for p in report.problems), (
        'The read noise should be checked against the nominal value'
    )
    assert 'REDUCTION NEEDS ATTENTION' in caplog.text, 'The verdict should be logged'
    assert 'Next: pipnick_astrometry' in caplog.text, 'The next step should be given'
    assert '[6/6]' in caplog.text, 'The progress should be logged'

    # The measurements of up-to-date products are read from their headers
    caplog.clear()
    report = StepReport('REDUCTION')
    reduce_all(night, rdxdir=rdxdir, report=report, qa=False)
    summary = Table.read(rdxdir / 'reduce' / 'summary.ecsv')
    assert list(summary['status']) == ['up to date'] * 6, 'Every product should be up to date'
    level = summary['bias_level'][summary['frametype'] == 'bias'][0]
    assert level == pytest.approx(300., abs=0.5), (
        'The measurements of up-to-date products should be reported'
    )
    assert len(report.mismatched) == 0, 'Products made by this version should not be reported'


def test_reduce_qa(night, tmp_path):
    rdxdir = tmp_path / 'rdx'
    qa_dir = rdxdir / 'reduce' / 'qa'
    reduce_all(night, rdxdir=rdxdir)
    built = {f.name: f.stat().st_mtime_ns for f in sorted(qa_dir.glob('*.png'))}
    products = sorted(f.stem for f in rdxdir.glob('*.fits'))
    assert sorted(f.removesuffix('.png') for f in built) == products, (
        'By default, a QA figure should be written for each product'
    )
    assert all(matplotlib.image.imread(f).ndim == 3 for f in sorted(qa_dir.glob('*.png'))), (
        'The QA figures should be readable PNG images'
    )

    # Up-to-date products only get a figure if they lack one
    (qa_dir / 'Bias.png').unlink()
    reduce_all(night, rdxdir=rdxdir)
    rerun = {f.name: f.stat().st_mtime_ns for f in sorted(qa_dir.glob('*.png'))}
    assert 'Bias.png' in rerun, 'A missing QA figure should be drawn'
    assert all(rerun[k] == built[k] for k in rerun if k != 'Bias.png'), (
        'The QA figures of up-to-date products should not be redrawn'
    )

    # Rebuilt products get a new figure
    reduce_all(night, rdxdir=rdxdir, overwrite=True)
    rebuilt = {f.name: f.stat().st_mtime_ns for f in sorted(qa_dir.glob('*.png'))}
    assert all(rebuilt[k] != rerun[k] for k in rerun), (
        'The QA figures of rebuilt products should be redrawn'
    )

    (qa_dir / 'Bias.png').unlink()
    reduce_all(night, rdxdir=rdxdir, qa=False)
    assert not (qa_dir / 'Bias.png').is_file(), 'No QA figure should be drawn if qa is False'


def test_reduce_version_mismatch(night, tmp_path, caplog):
    rdxdir = tmp_path / 'rdx'
    reduce_all(night, rdxdir=rdxdir, qa=False)
    with fits.open(rdxdir / 'DomeFlat_R.fits', mode='update') as hdu:
        hdu[0].header['PNVERS'] = '0.0.1'
    with fits.open(rdxdir / 'DomeFlat_V.fits', mode='update') as hdu:
        del hdu[0].header['PNVERS']
    report = StepReport('REDUCTION')
    reduce_all(night, rdxdir=rdxdir, report=report, qa=False)
    assert report.mismatched == {'DomeFlat_R.fits': '0.0.1', 'DomeFlat_V.fits': None}, (
        'Products made by another version, or with no version, should be reported'
    )
    assert 'made with a different version of pipnick' in caplog.text, (
        'The version mismatch should be reported in a warning'
    )
    assert report.status != 'FAILED', 'A version mismatch should not fail the step'
    assert fits.getheader(rdxdir / 'DomeFlat_R.fits')['PNVERS'] == '0.0.1', (
        'A version mismatch alone should not rebuild a product'
    )


def test_reduce_failure(night, tmp_path):
    rdxdir = tmp_path / 'rdx'
    # A corrupt dark frame
    (tmp_path / 'raw' / 'n5009.fits').write_bytes(b'not a FITS file')
    report = StepReport('REDUCTION')
    science = reduce_all(night, rdxdir=rdxdir, report=report, qa=False)
    summary = Table.read(rdxdir / 'reduce' / 'summary.ecsv')
    assert summary['status'][summary['file'] == 'Dark_120s.fits'][0] == 'failed', (
        'The failed product should be recorded in the summary'
    )
    assert science[0].is_file(), 'The other products should still be reduced'
    assert report.status == 'FAILED' and report.exit_status == 1, (
        'A failed product should fail the step, with a non-zero exit status'
    )


def test_reduction_status(night, tmp_path):
    rdxdir = tmp_path / 'rdx'
    reduce_all(night, rdxdir=rdxdir, qa=False)
    image = DataSet.from_file(night, rdxdir=rdxdir).science_images()[0]
    assert image.reduction_status() == ('ok', '-'), 'The reduced image should be up to date'
    set_column(night, ['n5011.fits', 'n5012.fits'], 'airmass', 1.5)
    image = DataSet.from_file(night, rdxdir=rdxdir).science_images()[0]
    assert image.reduction_status()[0] == 'stale', (
        'A change in the metadata table values should make the image stale'
    )
