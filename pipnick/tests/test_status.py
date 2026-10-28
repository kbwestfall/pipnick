"""
Tests for the status of the processing of a data set.
"""

from astropy.io import fits

from pipnick import metadata
from pipnick.pipelines.reduction import reduce_all
from pipnick.pipelines.status import night_status
from pipnick.tests.test_reduction import comb_id, set_column, set_comb_id, write_night


def build_night(tmp_path):
    """
    Build the per-setup file for the synthetic night of the reduction tests.
    """
    rawdir = tmp_path / 'raw'
    write_night(rawdir)
    files, _ = metadata.build_metadata(rawdir=rawdir, filename=tmp_path / 'n5.ecsv', split=True)
    table = files[0]
    set_comb_id(table, ['n5012.fits'], comb_id(table, 'n5011.fits'))
    return table


def test_status(tmp_path):
    table = build_night(tmp_path)
    rdxdir = tmp_path / 'rdx'
    before = table.read_text()
    tbl, next_step = night_status(table, rdxdir=rdxdir)
    assert table.read_text() == before, 'The status should never change the metadata file'
    assert list(tbl['reduced']) == ['missing'], 'An unreduced image should be missing'
    assert next_step.startswith('pipnick_reduce'), 'The next step should be the reduction'

    reduce_all(table, rdxdir=rdxdir)
    tbl, next_step = night_status(table, rdxdir=rdxdir)
    row = tbl[0]
    assert row['reduced'] == 'ok' and row['wcs'] == 'missing' and row['phot'] == 'missing', (
        'A reduced image without a solution or source table should be reported as such'
    )
    assert row['calibrated'] == 'n/a' and row['mosaic'] == 'n/a', (
        'Steps that do not apply should be n/a'
    )
    assert row['attention'] == '', 'An image made by this version needs no attention'
    assert next_step.startswith('pipnick_astrometry'), 'The next step should be the astrometry'

    # A change in the table values makes the image stale
    set_column(table, ['n5011.fits', 'n5012.fits'], 'airmass', 1.5)
    tbl, next_step = night_status(table, rdxdir=rdxdir)
    assert tbl['reduced'][0] == 'stale', 'The image should be stale'
    assert next_step.startswith('pipnick_reduce'), 'A stale image should be reduced again'

    # A reduction by another version needs attention
    path = rdxdir / str(tbl['image'][0])
    with fits.open(path, mode='update') as hdu:
        hdu[0].header['PNVERS'] = '0.0.1'
    tbl, _ = night_status(table, rdxdir=rdxdir)
    assert 'version' in tbl['attention'][0], 'A version mismatch should need attention'
