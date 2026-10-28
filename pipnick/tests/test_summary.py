"""
Tests for the summary tables of the steps.
"""

import pytest
from astropy.table import Table

from pipnick.pkg.summary import Summary
from pipnick.pkg.version import version


def test_summary(tmp_path):
    summary = Summary(
        'Test', ['image', 'status', 'nsrc', 'fwhm'], dtypes={'nsrc': int, 'fwhm': float},
        formats={'fwhm': '.2f'}, logged=['image', 'status']
    )
    row = summary.add(image='a.fits', status='failed')
    row['status'] = 'measured'
    summary.add(image='b.fits', status='measured', nsrc=10, fwhm=1.234)
    with pytest.raises(KeyError):
        summary.add(image='c.fits', bad=1)
    assert len(summary) == 2, 'The summary should have one row per item'

    tbl = summary.write(tmp_path / 'step')
    read = Table.read(tmp_path / 'step' / 'summary.ecsv')
    assert list(read['status']) == ['measured', 'measured'], (
        'Rows should be updatable until the summary is written'
    )
    assert bool(read['nsrc'].mask[0]) and not bool(read['nsrc'].mask[1]), (
        'Missing values should be masked'
    )
    assert read['nsrc'].dtype.kind == 'i' and read['fwhm'].dtype.kind == 'f', (
        'The columns should have their declared types'
    )
    assert read.meta['pipnick_version'] == version, 'The summary should record the version'
    assert tbl['fwhm'].info.format == '.2f', 'The column formats should be set'

    empty = Summary('Test', ['image'])
    assert empty.write(tmp_path / 'empty') is None, 'An empty summary should not be written'
    assert not (tmp_path / 'empty').exists(), 'An empty summary should not create a directory'
