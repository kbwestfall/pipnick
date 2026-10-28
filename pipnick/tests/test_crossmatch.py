"""
Tests for the cross-matching of celestial coordinates.
"""

import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from pipnick.pkg.exceptions import PipnickParameterError
from pipnick.utils.crossmatch import link_sources, match_sky

# 1 arcsec in degrees
ARCSEC = 1. / 3600.


def test_mutual():
    c1 = SkyCoord([10., 10.001, 10.002], [0., 0., 0.], unit='deg')
    c2 = SkyCoord([10.0000001, 10.00205, 11.], [0., 0., 0.], unit='deg')
    i1, i2, sep = match_sky(c1, c2, 1.)
    assert list(i1) == [0, 2] and list(i2) == [0, 1], (
        'Positions should be matched one-to-one within the radius'
    )
    assert np.all(sep < 1.), 'The separations should be within the radius'


def test_unique_second():
    # Two positions near the same position in the second list
    c1 = SkyCoord([10., 10. + 0.5 * ARCSEC, 10. + 5. * ARCSEC], [0., 0., 0.], unit='deg')
    c2 = SkyCoord([10. + 0.4 * ARCSEC, 10. + 5. * ARCSEC], [0., 0.], unit='deg')
    i1, i2, _ = match_sky(c1, c2, 1., unique='second')
    assert list(i1) == [1, 2] and list(i2) == [0, 1], (
        'Each position in the second list should be kept only in its nearest match'
    )
    i1, i2, _ = match_sky(c1, c2, 1., unique='none')
    assert list(i1) == [0, 1, 2] and list(i2) == [0, 0, 1], (
        'Without uniqueness, several positions can match the same position'
    )


def test_empty_and_options():
    c = SkyCoord([10.], [0.], unit='deg')
    i1, i2, sep = match_sky(SkyCoord([], [], unit='deg'), c, 1.)
    assert len(i1) == 0 and len(i2) == 0 and len(sep) == 0, (
        'Matching an empty list should return no matches'
    )
    with pytest.raises(PipnickParameterError):
        match_sky(c, c, 1., unique='first')


def test_link_sources():
    ra = [
        np.ma.MaskedArray([10., 10. + 10. * ARCSEC]),
        # Matches source 0, a new source, and a masked position
        np.ma.MaskedArray([10. + 0.3 * ARCSEC, 10. + 30. * ARCSEC, 10.], mask=[0, 0, 1]),
        # No coordinates
        np.ma.masked_all(2),
        # Matches sources 1 and 2
        np.ma.MaskedArray([10. + 10.2 * ARCSEC, 10. + 30.1 * ARCSEC])
    ]
    dec = [np.ma.zeros(len(r)) for r in ra]
    ids = link_sources(ra, dec, 1.)
    expected = [[0, 1], [0, 3, 4], [5, 6], [1, 3]]
    assert [list(i) for i in ids] == expected, (
        'Matched positions should take the identifier of their source; unmatched and masked '
        'positions should get new identifiers'
    )
