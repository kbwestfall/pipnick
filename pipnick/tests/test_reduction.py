"""
Tests for the reduction pipeline, using synthetic FITS files.
"""

import pytest

from pipnick import metadata
from pipnick.pipelines.reduction import reduce_all
from pipnick.pkg.exceptions import PipnickDataModelError
from pipnick.tests.test_metadata import write_andor


def test_one_setup(tmp_path):
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    write_andor(rawdir / 'n1001.fits')
    write_andor(rawdir / 'n1002.fits', AMPCONF='0.05MHz')
    table = tmp_path / 'raw.ecsv'
    metadata.build_metadata(rawdir=rawdir, filename=table)

    # Data from more than one setup cannot be reduced together
    with pytest.raises(PipnickDataModelError, match='more than one setup'):
        reduce_all(table=table, rdxdir=tmp_path / 'rdx')
