"""
Tests for the source catalog of a reduced image.
"""

import numpy as np
from astropy.io import fits
from astropy.table import MaskedColumn, Table
from astropy.wcs import WCS

from pipnick.datamodel.catalog import SourceCatalog


def catalog():
    """A small source table, with celestial coordinates."""
    tbl = Table()
    tbl['x'] = [10., 12., 50.]
    tbl['y'] = [10., 11., 50.]
    tbl['ra'] = MaskedColumn([150., 150.001, 150.01], mask=[False, False, True])
    tbl['dec'] = MaskedColumn([2., 2., 2.], mask=[False, False, True])
    tbl['flags'] = [0, 0, 1]
    tbl['flux_psf'] = MaskedColumn([1000., 50., 300.])
    tbl['flux_aper'] = MaskedColumn([1000., 50., 300.], mask=[False, True, False])
    tbl['flux_aper_err'] = MaskedColumn([10., 5., 3.])
    tbl['mag_inst'] = MaskedColumn([-7.5, -4.2, -6.2])
    tbl['mag_inst_err'] = MaskedColumn([0.01, 0.1, 0.02])
    tbl['iter_detected'] = [1, 2, 1]
    tbl.meta = {
        'exptime': 10., 'wcs_method': 'gaia', 'aper_radius_pix': 5.,
        'psf': {'fittype': 'circ', 'par': [0., 0., 1., 2., 3., 0.], 'fwhm_pix': 2.}
    }
    return SourceCatalog(tbl)


def test_access_and_selection():
    cat = catalog()
    assert len(cat) == 3 and list(cat['x']) == [10., 12., 50.], (
        'The columns should be accessed as those of the table'
    )
    assert cat.has_coordinates, 'The catalog should have celestial coordinates'
    assert list(cat.usable()) == [True, False, False], (
        'Usable sources should be unflagged, with an aperture flux and coordinates'
    )
    assert cat.psf.fittype == 'circ', 'The PSF model should come from the metadata'


def test_instrumental_mag():
    cat = catalog()
    mag, err = cat.instrumental_mag('psf')
    assert np.allclose(mag, [-7.5, -4.2, -6.2]), 'PSF magnitudes should be mag_inst'
    mag, err = cat.instrumental_mag('aper')
    assert np.isclose(mag[0], -2.5 * np.log10(100.)), (
        'Aperture magnitudes should be the aperture count rate'
    )
    assert bool(np.ma.getmaskarray(mag)[1]), 'Masked aperture fluxes should be masked'


def test_neighbors():
    cat = catalog()
    assert cat.neighbors(0) == (1, 0), (
        'A later detection within the residual radius should be a residual detection'
    )


def test_coordinates(tmp_path):
    cat = catalog()
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.crval = [150., 2.]
    wcs.wcs.crpix = [1., 1.]
    wcs.wcs.cd = [[-1e-4, 0.], [0., 1e-4]]
    hdr = wcs.to_header()
    hdr['PNASTMTH'] = 'gaia'
    cat.update_coordinates(fits.Header(hdr), 'abc')
    assert not np.any(np.ma.getmaskarray(cat['ra'])), 'Every source should have coordinates'
    assert cat.meta['wcs_hash'] == 'abc', 'The solution hash should be recorded'
    cat.write(tmp_path / 'phot.ecsv')
    read = SourceCatalog.read(tmp_path / 'phot.ecsv')
    assert np.allclose(read['ra'], cat['ra']) and read.path == tmp_path / 'phot.ecsv', (
        'The catalog should be read back'
    )
