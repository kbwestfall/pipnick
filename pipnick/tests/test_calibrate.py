"""
Tests for the photometric calibration, using synthetic photometry tables of
Landolt standards with known calibration coefficients.
"""

import logging
import tomllib

import numpy as np
import pytest
from astropy import units
from astropy.io import fits
from astropy.table import MaskedColumn, Table
from astropy.time import Time

from pipnick.metadata.build import build_metadata
from pipnick.astrometry.services import propagate_proper_motion
from pipnick.datamodel.calibration import Calibration, toml_value
from pipnick.datamodel.catalog import SourceCatalog
from pipnick.datamodel.dataset import DataSet
from pipnick.photometry import extinction, landolt, nightfit
from pipnick.photometry.calibrate import (
    CALIBRATION_TERMS, check_terms, design_matrix, fit_calibration, iterate_colors,
    residual_drift, zero_point
)
from pipnick.pipelines import calibrate
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickParameterError
from pipnick.pkg.provenance import image_hash
from pipnick.pkg.version import version
from pipnick.tests.test_metadata import write_andor

# The standards: well-observed stars of Landolt field SA 110
STARS = ['110 229', '110 230', '110 232', '110 233', '110 340', '110 477', '110 355', '110 361']
# The injected calibration of each band: zero point, extinction, color term
TRUE = {'B': (-24.1, 0.30, 0.05), 'V': (-24.5, 0.17, -0.03), 'R': (-24.6, 0.12, 0.02)}
NOISE = 0.004
EXPTIME = 30.
# The FWHM of the PSF (pixels)
FWHM = 4.


def read_calibration(path):
    """
    Read a calibration file, and return its content and whether it has been
    edited.
    """
    calibration = Calibration.read(path)
    return calibration.content, calibration.edited


def standards():
    """
    Return the catalog entries of the standards.
    """
    catalog = landolt.load_landolt()
    return catalog[[list(catalog['name']).index(s) for s in STARS]]


def instrumental(mag, color, airmass, band, rng):
    """
    Return instrumental magnitudes for the injected calibration.
    """
    k0, k1, kc = TRUE[band]
    return mag + k0 + k1 * airmass + kc * color + rng.normal(0., NOISE, np.size(mag))


def phot_table(ra, dec, mag, meta):
    """
    Construct a photometry source table; see
    :func:`pipnick.pipelines.photometry.photometry_all`.
    """
    n = len(ra)
    nomask = np.zeros(n, dtype=bool)
    flux = EXPTIME * 10**(-0.4 * np.asarray(mag))
    tbl = Table()
    tbl['id'] = np.arange(1, n + 1)
    tbl['group_id'] = np.arange(1, n + 1)
    tbl['group_size'] = np.ones(n, dtype=int)
    tbl['iter_detected'] = np.ones(n, dtype=int)
    tbl['x'] = MaskedColumn(np.linspace(10., 200., n), mask=nomask)
    tbl['y'] = MaskedColumn(np.linspace(20., 210., n), mask=nomask)
    tbl['ra'] = MaskedColumn(np.asarray(ra, dtype=float), mask=nomask, unit='deg')
    tbl['dec'] = MaskedColumn(np.asarray(dec, dtype=float), mask=nomask, unit='deg')
    tbl['flux_psf'] = MaskedColumn(flux, mask=nomask)
    tbl['flux_aper'] = MaskedColumn(flux, mask=nomask)
    tbl['flux_aper_err'] = MaskedColumn(flux * NOISE / 1.0857, mask=nomask)
    tbl['mag_inst'] = MaskedColumn(np.asarray(mag, dtype=float), mask=nomask)
    tbl['mag_inst_err'] = MaskedColumn(np.full(n, NOISE), mask=nomask)
    tbl['flags'] = np.zeros(n, dtype=int)
    tbl.meta = dict(meta, psf={'fwhm_pix': FWHM}, aper_radius_pix=2 * FWHM)
    return tbl


def write_night(tmp_path, design, seed=1, science=False, outlier=None):
    """
    Write the metadata table, reduced images, and photometry tables of a
    night of observations of the standards.

    Parameters
    ----------
    design : list
        The filter, airmass, and UT hour of each standard observation.
    science : bool
        Also observe a science field (``M11``) in B and V, with six stars.
    outlier : str
        A star whose instrumental magnitudes are offset by 0.3 mag.

    Returns
    -------
    table, rdxdir : :class:`pathlib.Path`
        The per-setup table and the reduction directory.
    truth : dict
        The true magnitudes of the science stars.
    """
    rng = np.random.default_rng(seed)
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    frames = [('SA110', f, a, h) for f, a, h in design]
    if science:
        frames += [('M11', 'B', 1.3, 6.), ('M11', 'V', 1.32, 6.2), ('M11', 'V', 1.34, 6.4)]
    for i, (obj, filt, airmass, hour) in enumerate(frames):
        date = Time('2026-08-20T00:00:00', scale='utc') + hour * units.hour
        write_andor(
            rawdir / f'n9{i:03d}.fits', OBJECT=f'{obj} {filt}', FILTER=filt, AIRMASS=airmass,
            XPOSURE=EXPTIME, **{'DATE-BEG': date.isot}
        )
    files, _ = build_metadata(rawdir=rawdir, filename=tmp_path / 'n9.ecsv', split=True)
    table = files[0]
    rdxdir = tmp_path / 'rdx'
    (rdxdir / 'photometry').mkdir(parents=True)
    dataset = DataSet.from_file(table, rdxdir=rdxdir)
    tbl, tag = dataset.table, dataset.tag

    stars = standards()
    sci = {
        'ra': 282.77 + np.arange(6) * 0.003, 'dec': np.full(6, -6.27),
        'V': np.linspace(12., 15., 6), 'B-V': np.linspace(0.2, 1.4, 6)
    }
    truth = {'V': sci['V'], 'B': sci['V'] + sci['B-V'], 'B-V': sci['B-V']}
    for product in dataset.products():
        i = product.indx[0]
        obj, filt, airmass, _ = frames[int(str(tbl['file'][i])[2:5])]
        mjd = float(tbl['mjd'][i])
        hdr = fits.Header(
            {
                'PNSETUP': tag, 'PNRAWHSH': f'raw{product.comb_id}', 'PNTABH': 'table',
                'PNASTMTH': 'gaia', 'PNASTWH': f'wcs{product.comb_id}', 'AIRMASS': airmass
            }
        )
        fits.PrimaryHDU(data=np.zeros((8, 8)), header=hdr).writeto(product.file)
        meta = {
            'image': product.file.name, 'setup': tag, 'comb_id': product.comb_id,
            'object': obj, 'filter': filt, 'exptime': EXPTIME, 'airmass': airmass, 'mjd': mjd,
            'settings_hash': 'settings', 'image_hash': image_hash(hdr), 'wcs_method': 'gaia',
            'wcs_hash': f'wcs{product.comb_id}'
        }
        if obj == 'SA110':
            epoch = Time(mjd, format='mjd').jyear
            coo = propagate_proper_motion(stars, epoch)
            # Filters other than UBVRI are given the V calibration
            band = filt if filt in TRUE else 'V'
            color, _ = landolt.catalog_color(stars, calibrate.default_colors()[band])
            mag = instrumental(np.asarray(stars[band]), np.asarray(color), airmass, band, rng)
            if outlier is not None:
                mag[STARS.index(outlier)] += 0.3
            # Add sources that are not standards
            ra = np.append(coo.ra.deg, coo.ra.deg[0] + 0.01)
            dec = np.append(coo.dec.deg, coo.dec.deg[0])
            mag = np.append(mag, -10.)
        else:
            ra, dec = sci['ra'], sci['dec']
            mag = instrumental(truth[filt], sci['B-V'], airmass, filt, rng)
        phot_table(ra, dec, mag, meta).write(
            rdxdir / 'photometry' / f'{product.file.stem}_phot.ecsv', format='ascii.ecsv'
        )
    return table, rdxdir, truth


def wide_design():
    """
    Return a design with standards over a wide range of airmass in B, V, R.
    """
    return [(f, a, h) for h, a in [(4., 1.1), (5., 1.5), (7., 2.0)] for f in 'BVR']


def narrow_design():
    """
    Return a design with standards over a narrow range of airmass in B, V.
    """
    return [(f, a, h) for h, a in [(4., 1.18), (5., 1.25)] for f in 'BV']


def calibration_file(rdxdir):
    """
    Return the calibration file of the synthetic night.
    """
    return sorted((rdxdir / 'calibrate').glob('*.toml'))[0]


def test_landolt_catalog():
    catalog = landolt.load_landolt()
    assert len(catalog) == 942, 'The packaged catalog should merge the 2009 and 2013 catalogs'
    assert np.all(catalog['ref_epoch'] == 2000.), 'The positions should be at epoch 2000'
    good = np.logical_not(
        np.ma.getmaskarray(np.ma.asarray(catalog['B']))
        | np.ma.getmaskarray(np.ma.asarray(catalog['B-V']))
    )
    assert np.allclose(
        np.asarray(catalog['B'] - catalog['V'])[good], np.asarray(catalog['B-V'])[good],
        atol=1e-3
    ), 'The derived B magnitudes should be consistent with V and B-V'
    good = np.logical_not(np.ma.getmaskarray(np.ma.asarray(catalog['I'])))
    assert np.allclose(
        np.asarray(catalog['V'] - catalog['I'])[good], np.asarray(catalog['V-I'])[good],
        atol=1e-3
    ), 'The derived I magnitudes should be consistent with V and V-I'
    color, _ = landolt.catalog_color(catalog[:5], 'B-R')
    assert np.allclose(color, catalog['B'][:5] - catalog['R'][:5]), (
        'An unpublished color should be the difference of the derived magnitudes'
    )


def test_parse_color():
    assert landolt.parse_color('B-V') == ('B', 'V'), 'A color should give its two bands'
    for bad in ['BV', 'B-B', 'B-g', 'B-V-R']:
        with pytest.raises(PipnickParameterError):
            landolt.parse_color(bad)


def test_match_landolt():
    stars = standards()[:3]
    coo = propagate_proper_motion(stars, 2026.)
    # Two sources near the first star, one 1.5 arcsec away; one far source
    ra = np.array([coo.ra.deg[0], coo.ra.deg[0] + 1.5 / 3600, coo.ra.deg[1], 10.])
    dec = np.array([coo.dec.deg[0], coo.dec.deg[0], coo.dec.deg[1], 10.])
    isrc, icat, sep = landolt.match_landolt(stars, ra, dec, 2026., radius=2.)
    assert list(isrc) == [0, 2] and list(icat) == [0, 1], (
        'Each standard should be matched to its nearest source, and the far source not at all'
    )
    assert np.all(sep < 0.01), 'The propagated positions should match exactly'
    _, _, sep = landolt.match_landolt(stars, ra, dec, 2000., radius=2.)
    assert sep[0] > 0.01, 'Ignoring the proper motion should increase the separation'


def test_band_extinction(caplog):
    wave, ext = extinction.load_extinction(extinction.extinction_path('mthamextinct.dat'))
    expected = {'U': 0.604, 'B': 0.300, 'V': 0.178, 'R': 0.120, 'I': 0.075}
    with caplog.at_level(logging.WARNING):
        values = {b: extinction.band_extinction(wave, ext, b) for b in expected}
    for band, k in expected.items():
        assert values[band][0] == pytest.approx(k, abs=2e-3), (
            f'The nominal {band} extinction should match the folded Mt. Hamilton curve'
        )
    assert values['U'][1] < 1 and values['V'][1] == 1., (
        'The curve should only partially cover the U band'
    )
    assert 'covers only' in caplog.text, 'The partial coverage of U should be reported'


def test_fit_calibration():
    rng = np.random.default_rng(2)
    n = 40
    airmass = rng.uniform(1., 2., n)
    color = rng.uniform(0., 1.5, n)
    mag_c = rng.uniform(11., 15., n)
    k = np.array([-24., 0.2, 0.05, 0.01, -0.02])
    x = design_matrix(CALIBRATION_TERMS, airmass, color)
    mag_i = mag_c + np.dot(k, x) + rng.normal(0., 0.002, n)
    err = np.full(n, 0.002)

    fit = fit_calibration(
        mag_c, color, mag_i, airmass, list(CALIBRATION_TERMS), err=err
    )
    assert np.all(np.abs(fit.coeff - k) < 4 * fit.errors + 1e-6), (
        'All of the coefficients should be recovered within their errors'
    )
    # No mask given: all measurements are used (the old code built a mask of
    # the wrong shape)
    assert fit.gpm.shape == (n,) and np.all(fit.gpm), 'Without a mask, all should be used'

    fit = fit_calibration(
        mag_c, color, mag_i, airmass, ['zp', 'airmass', 'color', 'color2', 'airmass_color'],
        err=err, fixed={'airmass': 0.2}
    )
    assert fit.coeff[1] == 0.2 and fit.errors[1] == 0., (
        'A fixed term should keep its value and have no error'
    )
    assert np.all(fit.cov[1] == 0.), 'A fixed term should not contribute to the covariance'

    with pytest.raises(PipnickParameterError):
        fit_calibration(
            mag_c[:3], color[:3], mag_i[:3], airmass[:3], ['zp', 'color', 'airmass']
        )
    with pytest.raises(PipnickParameterError):
        check_terms(['airmass'])


def test_zero_point():
    terms = ['zp', 'airmass', 'color', 'airmass_color']
    coeff = np.array([-24., 0.2, 0.05, -0.02])
    cov = np.diag([1e-4, 4e-4, 1e-4, 0.])
    zp = zero_point(terms, coeff, cov, 1.5)
    assert zp['zp'] == pytest.approx(24. - 0.3), 'The zero point should include the extinction'
    assert zp['c1'] == pytest.approx(-(0.05 - 0.03)), (
        'The color coefficient should include the airmass-color term'
    )
    assert zp['zp_err'] == pytest.approx(np.sqrt(1e-4 + 1.5**2 * 4e-4)), (
        'The zero-point error should be propagated from the covariance'
    )


def test_residual_drift():
    rng = np.random.default_rng(3)
    time = np.linspace(0., 8., 30)
    err = np.full(30, 0.005)
    slope, slope_err = residual_drift(time, 0.01 * time + rng.normal(0., 0.005, 30), err)
    assert abs(slope - 0.01) < 3 * slope_err and slope_err < 0.002, (
        'An injected drift should be recovered'
    )
    slope, slope_err = residual_drift(time, rng.normal(0., 0.005, 30), err)
    assert abs(slope) < 3 * slope_err, 'No drift should be found without one'
    assert residual_drift(time[:2], time[:2], err[:2]) is None, (
        'A drift should not be measured from fewer than three measurements'
    )


def test_iterate_colors():
    terms = ['zp', 'airmass', 'color']
    cals = {
        'B': {
            'terms': terms, 'coeff': np.array([-24.1, 0.3, 0.05]), 'cov': np.zeros((3, 3)),
            'color': ('B', 'V')
        },
        'V': {
            'terms': terms, 'coeff': np.array([-24.5, 0.17, -0.03]), 'cov': np.zeros((3, 3)),
            'color': ('B', 'V')
        },
        'R': {
            'terms': terms, 'coeff': np.array([-24.6, 0.12, 0.02]), 'cov': np.zeros((3, 3)),
            'color': ('V', 'R')
        }
    }
    # Source 0: B, V (twice) and R; source 1: B and V only; source 2: V only
    src = np.array([0, 0, 0, 0, 1, 1, 2])
    band = np.array(['B', 'V', 'V', 'R', 'B', 'V', 'V'])
    airmass = np.array([1.2, 1.3, 1.5, 1.4, 1.2, 1.3, 1.3])
    true = {
        (0, 'B'): 13.0, (0, 'V'): 12.2, (0, 'R'): 11.7, (1, 'B'): 14.0, (1, 'V'): 13.5,
        (2, 'V'): 12.0
    }
    true_color = {0: {'B': 0.8, 'V': 0.8, 'R': 0.5}, 1: {'B': 0.5, 'V': 0.5}, 2: {'V': 0.4}}
    mag_i = np.array(
        [
            true[(s, b)] + np.dot(cals[b]['coeff'], [1., a, true_color[s][b]])
            for s, b, a in zip(src, band, airmass)
        ]
    )
    mag, mag_err, color, assumed = iterate_colors(
        np.ma.MaskedArray(mag_i), np.ma.MaskedArray(np.full(7, 0.01)), airmass, src, band, cals,
        default_color={('B', 'V'): 0.4}
    )
    expected = np.array([true[(s, b)] for s, b in zip(src, band)])
    assert np.allclose(mag[:6], expected[:6], atol=1e-3), (
        'The magnitudes of sources with measured colors should be recovered'
    )
    assert not np.any(assumed[:6]), (
        'A color should be measured when both of its bands are, even if another band is missing'
    )
    assert assumed[6] and color[6] == pytest.approx(0.4), (
        'A source seen in one band should use the default color'
    )
    assert mag[6] == pytest.approx(12.0, abs=1e-3), (
        'The default color equal to the true color should give the true magnitude'
    )


def test_calibrate_all(tmp_path, caplog):
    table, rdxdir, _ = write_night(tmp_path, wide_design())
    with caplog.at_level(logging.INFO):
        outputs = calibrate.calibrate_all(table, rdxdir=rdxdir, qa=True)
    assert len(outputs) == 1, 'There should be one calibration for the night'
    content, edited = read_calibration(outputs[0])
    assert not edited, 'A new calibration file should not be reported as edited'
    assert content['calibration']['night'] == '2026-08-20', 'The night should be the UT date'
    for band, (k0, k1, kc) in TRUE.items():
        entry = content['band'][band]
        assert entry['status'] == 'ok', f'Band {band} should be calibrated'
        assert entry['source'] == ['fit', 'fit', 'fit'], 'All terms should be fit'
        k = np.array(entry['coefficients'])
        assert np.all(np.abs(k - [k0, k1, kc]) < 4 * np.array(entry['errors']) + 0.003), (
            f'The coefficients of band {band} should be recovered'
        )
        assert entry['n_stars'] == len(STARS), 'All of the standards should be used'
        assert entry['rms'] < 3 * NOISE, 'The residuals should be consistent with the noise'
        assert entry['n_loo'] == len(STARS), 'Every star should have a leave-one-out residual'
    resid = Table.read(sorted((rdxdir / 'calibrate').glob('*_residuals.ecsv'))[0])
    assert len(resid) == 9 * len(STARS), 'Every standard measurement should be in the table'
    assert np.all(resid['used']), 'All of the standard measurements should be used'
    assert resid.meta['extinction'] == {'B': 'fit', 'V': 'fit', 'R': 'fit'}, (
        'The residual table should record how the extinction was set'
    )
    assert (rdxdir / 'calibrate' / 'summary.ecsv').is_file(), 'The summary should be written'
    assert len(sorted((rdxdir / 'calibrate' / 'qa').glob('*.png'))) == 1, (
        'The QA figure should be written'
    )
    assert 'EXTINCTION NOT MEASURED' not in caplog.text, 'There should be no fallback'

    # Rerun: up to date
    caplog.clear()
    text = outputs[0].read_text()
    with caplog.at_level(logging.INFO):
        calibrate.calibrate_all(table, rdxdir=rdxdir)
    assert 'is up to date' in caplog.text and outputs[0].read_text() == text, (
        'A rerun with unchanged inputs should not refit'
    )

    # Edit the calibration: it is kept
    caplog.clear()
    outputs[0].write_text(text.replace('[band.V]\n', '[band.V]\n# edited\n'))
    edited_text = outputs[0].read_text().replace(
        f'coefficients = {toml_value(content["band"]["V"]["coefficients"])}',
        'coefficients = [-24.0, 0.17, -0.03]'
    )
    outputs[0].write_text(edited_text)
    with caplog.at_level(logging.INFO):
        calibrate.calibrate_all(table, rdxdir=rdxdir)
    assert 'has been edited; keeping it' in caplog.text, 'The edit should be detected'
    assert outputs[0].read_text() == edited_text, 'An edited file should be kept'

    # Changed inputs with an edited file: a warning
    caplog.clear()
    with caplog.at_level(logging.INFO):
        calibrate.calibrate_all(table, rdxdir=rdxdir, sys_err=0.02)
    assert 'may hide changes' in caplog.text, (
        'Changed inputs of an edited file should give a warning'
    )
    assert outputs[0].read_text() == edited_text, 'The edited file should still be kept'

    # Overwrite refits
    calibrate.calibrate_all(table, rdxdir=rdxdir, overwrite=True)
    content, edited = read_calibration(outputs[0])
    assert not edited and content['band']['V']['coefficients'][0] != -24.0, (
        '--overwrite should refit an edited file'
    )


def test_fallback(tmp_path, caplog):
    table, rdxdir, _ = write_night(tmp_path, narrow_design())
    with caplog.at_level(logging.WARNING):
        calibrate.calibrate_all(
            table, rdxdir=rdxdir, terms=['zp', 'airmass', 'color', 'airmass_color'],
            apply='image'
        )
    cfile = calibration_file(rdxdir)
    content, _ = read_calibration(cfile)
    assert content['calibration']['pipnick_version'] == version, (
        'The calibration file should record the pipnick version'
    )
    nominal = extinction.nominal_extinction(
        DataSet.from_file(table, rdxdir=rdxdir).camera, ['B', 'V']
    )
    for band in ['B', 'V']:
        entry = content['band'][band]
        assert entry['source'] == ['fit', 'fallback', 'fit', 'fallback'], (
            'A narrow airmass span should fall back to the nominal extinction and drop the '
            'airmass-color term'
        )
        assert entry['coefficients'][1] == pytest.approx(nominal[band]), (
            'The fallback should use the nominal extinction'
        )
        assert entry['coefficients'][3] == 0., 'The airmass-color term should be fixed at 0'
        assert entry['errors'][1] == 0. and entry['errors'][3] == 0., (
            'Fallback terms should have no error'
        )
        assert 'airmass span 0.07 < 0.2' == entry['fallback_reason'], (
            'The reason for the fallback should be recorded'
        )
    text = cfile.read_text()
    assert text.count('extinction NOT MEASURED') == 2, (
        'The calibration file should say that the extinction was not measured'
    )
    hdr = fits.getheader(sorted(rdxdir.glob('SA110_B_*.fits'))[0])
    assert hdr['PNCALEXT'] == 'fallback', 'The header should record the fallback'
    resid = Table.read(sorted((rdxdir / 'calibrate').glob('*_residuals.ecsv'))[0])
    assert resid.meta['extinction'] == {'B': 'fallback', 'V': 'fallback'}, (
        'The residual table should record the fallback'
    )
    summary = Table.read(rdxdir / 'calibrate' / 'summary.ecsv')
    assert list(summary['extinction']) == ['fallback', 'fallback'], (
        'The summary should record the fallback'
    )
    assert 'EXTINCTION NOT MEASURED: the nominal extinction was used' in caplog.text, (
        'The run should end with a banner listing the fallbacks'
    )

    # Without the fallback, every band fails
    with pytest.raises(PipnickDataModelError, match='No band could be calibrated'):
        calibrate.calibrate_all(table, rdxdir=rdxdir, fallback=False, overwrite=True)
    content, _ = read_calibration(cfile)
    assert content['band']['B']['status'] == 'failed', 'The bands should be recorded as failed'
    assert 'fallback is disabled' in content['band']['B']['reason'], (
        'The reason for the failure should be recorded'
    )


def test_sample_checks(tmp_path, caplog):
    design = [('B', 1.1, 4.), ('B', 1.6, 5.), ('V', 1.1, 4.), ('V', 1.6, 5.)]
    table, rdxdir, _ = write_night(tmp_path, design)
    calibrate.calibrate_all(table, rdxdir=rdxdir)
    content, _ = read_calibration(calibration_file(rdxdir))
    assert content['band']['B']['status'] == 'ok', 'Eight standards should be sufficient'

    with caplog.at_level(logging.WARNING):
        calibrate.calibrate_all(table, rdxdir=rdxdir, min_stars=10, overwrite=True, force=True)
    content, _ = read_calibration(calibration_file(rdxdir))
    assert content['band']['B']['status'] == 'ok' and 'fewer than 10' in caplog.text, (
        '--force should downgrade the minimum-stars error to a warning'
    )
    with pytest.raises(PipnickDataModelError):
        calibrate.calibrate_all(table, rdxdir=rdxdir, min_stars=10, overwrite=True)
    content, _ = read_calibration(calibration_file(rdxdir))
    assert 'fewer than 10' in content['band']['B']['reason'], (
        'Too few stars should fail the band'
    )
    with pytest.raises(PipnickDataModelError):
        calibrate.calibrate_all(table, rdxdir=rdxdir, min_color_span=5., overwrite=True)
    content, _ = read_calibration(calibration_file(rdxdir))
    assert 'color span' in content['band']['B']['reason'], (
        'A small color span should fail the band'
    )


def test_band_failure_isolated(tmp_path):
    # B has one image, with only three usable standards; V has two images
    design = [('B', 1.1, 4.), ('V', 1.1, 4.), ('V', 1.6, 5.)]
    table, rdxdir, _ = write_night(tmp_path, design)
    pfile = sorted((rdxdir / 'photometry').glob('SA110_B_*_phot.ecsv'))[0]
    phot = Table.read(pfile, format='ascii.ecsv')
    phot['flags'][3:len(STARS)] = 1
    phot.write(pfile, format='ascii.ecsv', overwrite=True)
    terms = ['zp', 'airmass', 'color', 'color2']
    calibrate.calibrate_all(
        table, rdxdir=rdxdir, terms=terms, min_airmass_span=0., min_color_span=0., force=True,
        fallback=False
    )
    content, _ = read_calibration(calibration_file(rdxdir))
    assert content['band']['V']['status'] == 'ok', 'V should be calibrated'
    assert content['band']['B']['status'] == 'failed', (
        'A band with too few measurements for its terms should fail, even with --force'
    )
    assert 'only 3 measurements for 4 fitted terms' in content['band']['B']['reason'], (
        'The reason should be the number of measurements'
    )


def test_degenerate_terms(tmp_path):
    # B has a single airmass, so its zero point and extinction are degenerate
    design = [('B', 1.1, 4.), ('V', 1.1, 4.), ('V', 1.6, 5.)]
    table, rdxdir, _ = write_night(tmp_path, design)
    calibrate.calibrate_all(table, rdxdir=rdxdir, min_airmass_span=0., fallback=False)
    content, _ = read_calibration(calibration_file(rdxdir))
    assert content['band']['V']['status'] == 'ok', 'V should be calibrated'
    assert 'degenerate' in content['band']['B']['reason'], (
        'Degenerate terms should fail the band instead of giving an invalid fit'
    )


def add_source(phot, i, dx, flux_ratio, iteration):
    """
    Add a source near source ``i`` of a photometry table.
    """
    phot.add_row(phot[i])
    phot['id'][-1] = len(phot)
    phot['x'][-1] += dx
    phot['flux_psf'][-1] *= flux_ratio
    phot['iter_detected'][-1] = iteration


def test_neighbors():
    meta = {'psf': {'fwhm_pix': FWHM}, 'aper_radius_pix': 2 * FWHM}
    phot = phot_table(np.zeros(1), np.zeros(1), np.full(1, -12.), meta)
    add_source(phot, 0, 0.5 * FWHM, 0.03, 2)
    assert SourceCatalog(phot).neighbors(0) == (1, 0), (
        'A faint later detection very close to a source should be a residual detection'
    )
    add_source(phot, 0, 1.2 * FWHM, 0.05, 1)
    assert SourceCatalog(phot).neighbors(0) == (1, 1), (
        'A first-iteration neighbor within the aperture should contaminate the source'
    )
    add_source(phot, 0, 1.9 * FWHM, 0.005, 1)
    add_source(phot, 0, 3 * FWHM, 1., 1)
    assert SourceCatalog(phot).neighbors(0) == (1, 1), (
        'Faint neighbors and neighbors beyond the aperture should not contaminate the source'
    )
    assert SourceCatalog(phot).neighbors(0, max_neighbor_flux=0.001) == (1, 2), (
        'A lower flux limit should make the faint neighbor contaminating'
    )


def test_selection(tmp_path):
    table, rdxdir, _ = write_night(tmp_path, wide_design())
    pfile = sorted((rdxdir / 'photometry').glob('SA110_V_*_phot.ecsv'))[0]
    phot = Table.read(pfile, format='ascii.ecsv')
    phot['flags'][0] = 4
    mask = np.zeros(len(phot), dtype=bool)
    mask[1] = True
    phot['flux_aper'] = MaskedColumn(phot['flux_aper'], mask=mask)
    phot['ra'][3] += 5. / 3600
    # A contaminating neighbor of star 2, and a residual detection near star 4
    add_source(phot, 2, 1.2 * FWHM, 0.05, 1)
    add_source(phot, 4, 0.5 * FWHM, 0.03, 2)
    phot['group_size'][4] = 2
    phot.write(pfile, format='ascii.ecsv', overwrite=True)
    calibrate.calibrate_all(table, rdxdir=rdxdir)
    resid = Table.read(sorted((rdxdir / 'calibrate').glob('*_residuals.ecsv'))[0])
    resid = resid[resid['image'] == pfile.name.replace('_phot.ecsv', '')]
    reasons = dict(zip(resid['star'], resid['reason']))
    assert reasons[STARS[0]] == 'flagged', 'A flagged source should not be used'
    assert reasons[STARS[1]] == 'saturated or otherwise masked', (
        'A source with masked aperture pixels should not be used'
    )
    assert reasons[STARS[2]] == 'contaminated', 'A contaminated source should not be used'
    assert STARS[3] not in reasons, 'A source beyond the match radius should not be matched'
    assert reasons[STARS[4]] == '-', (
        'A source fit together with a residual detection should still be used'
    )
    assert resid['n_residual'][resid['star'] == STARS[4]][0] == 1, (
        'The residual detection should be counted'
    )
    assert not np.any(resid['used'][resid['reason'] != '-']), (
        'Unusable measurements should not be used'
    )


def test_inputs(tmp_path, caplog):
    table, rdxdir, _ = write_night(tmp_path, wide_design())
    pfiles = sorted((rdxdir / 'photometry').glob('*_phot.ecsv'))
    # A missing table is skipped with a warning
    pfiles[0].unlink()
    with caplog.at_level(logging.WARNING):
        calibrate.calibrate_all(table, rdxdir=rdxdir)
    assert 'has no source table; run pipnick_photometry' in caplog.text, (
        'A missing source table should be reported'
    )
    # A stale table is an error
    phot = Table.read(pfiles[1], format='ascii.ecsv')
    good_hash = phot.meta['image_hash']
    phot.meta['image_hash'] = 'old'
    phot.write(pfiles[1], format='ascii.ecsv', overwrite=True)
    with pytest.raises(PipnickDataModelError, match='out of date'):
        calibrate.calibrate_all(table, rdxdir=rdxdir)
    # Tables without coordinates are an error
    for pfile in pfiles[1:]:
        phot = Table.read(pfile, format='ascii.ecsv')
        phot.meta['wcs_method'] = 'none'
        if pfile == pfiles[1]:
            phot.meta['image_hash'] = good_hash
        phot.write(pfile, format='ascii.ecsv', overwrite=True)
    with pytest.raises(PipnickDataModelError, match='No source table has celestial'):
        calibrate.calibrate_all(table, rdxdir=rdxdir)


def test_apply(tmp_path, caplog):
    table, rdxdir, truth = write_night(tmp_path, wide_design(), science=True)
    calibrate.calibrate_all(table, rdxdir=rdxdir, apply='both')

    # The header keywords reproduce the catalog magnitudes of the standards
    stars = standards()
    path = sorted(rdxdir.glob('SA110_V_*.fits'))[0]
    hdr = fits.getheader(path)
    phot = Table.read(rdxdir / 'photometry' / f'{path.stem}_phot.ecsv', format='ascii.ecsv')
    color = np.asarray(stars['B-V'])
    mag = np.asarray(phot['mag_inst'][:len(STARS)]) + hdr['PNCALZP'] + hdr['PNCALCT'] * color
    assert np.allclose(mag, np.asarray(stars['V']), atol=4 * NOISE), (
        'The header calibration should reproduce the catalog magnitudes'
    )
    assert hdr['PNCALEXT'] == 'fit' and hdr['PNCALCOL'] == 'B-V', (
        'The header should record the extinction source and the color'
    )
    assert hdr['PHOTZP'] == hdr['PNCALZP'], 'PHOTZP should equal PNCALZP'
    assert hdr['PNCALVER'] == version, 'The calibration should record the pipnick version'
    assert np.all(fits.getdata(path) == 0.), 'The pixels should not be changed'
    assert image_hash(hdr) == phot.meta['image_hash'], (
        'The photometry should still see the image as up to date'
    )

    # A rerun with an unchanged calibration leaves the header alone
    caplog.clear()
    with caplog.at_level(logging.INFO):
        calibrate.calibrate_all(table, rdxdir=rdxdir, apply='image')
    assert 'header calibration is up to date' in caplog.text, (
        'An unchanged calibration should not rewrite the headers'
    )

    # The calibrated source tables of the science field
    tables = {
        f.name: Table.read(f) for f in sorted((rdxdir / 'calibrate').glob('M11_*_cal.ecsv'))
    }
    assert len(tables) == 3, 'Each science image should have a calibrated table'
    for name, tbl in tables.items():
        band = tbl.meta['band']
        assert np.allclose(tbl['mag'], truth[band], atol=5 * NOISE), (
            f'The calibrated {band} magnitudes should be recovered'
        )
        assert not np.any(tbl['color_assumed']), 'The colors should be measured'
        assert np.allclose(tbl['color_used'], truth['B-V'], atol=5 * NOISE), (
            'The measured colors should be recovered'
        )
        assert len(set(tbl['source_id'])) == 6, 'The sources should keep their identifiers'
    ids = [list(tbl['source_id']) for tbl in tables.values()]
    assert ids[0] == ids[1] == ids[2], 'The same sources should be matched across the images'


def test_leave_one_out(tmp_path):
    # The outlier has a color in the middle of the range, so that it does not
    # drive the color term
    table, rdxdir, _ = write_night(tmp_path, wide_design(), outlier=STARS[6])
    calibrate.calibrate_all(table, rdxdir=rdxdir, rej=100.)
    resid = Table.read(sorted((rdxdir / 'calibrate').glob('*_residuals.ecsv'))[0])
    resid = resid[resid['band'] == 'V']
    out = resid['star'] == STARS[6]
    assert np.all(np.abs(resid['resid_loo'][out]) > 0.25), (
        'The outlier should stand out in its leave-one-out residuals'
    )
    assert np.all(np.abs(resid['resid_loo'][out]) > np.abs(resid['resid'][out])), (
        'The leave-one-out residuals of an outlier should exceed its residuals'
    )
    # The outlier shifts the zero point of the fits without each of the other
    # stars by about 0.3/8 mag
    assert np.min(np.abs(resid['resid_loo'][out])) > 3 * np.max(
        np.abs(resid['resid_loo'][np.logical_not(out)])
    ), 'The other stars should have much smaller leave-one-out residuals'


def test_leave_one_out_no_dof(tmp_path):
    loo = nightfit._leave_one_out
    meas = Table(
        {
            'star': ['a', 'b', 'c'], 'airmass': [1.1, 1.5, 2.], 'color': [0.1, 0.5, 1.],
            'mag_cat': [12., 13., 14.], 'mag_inst': [-12., -11., -10.], 'err': [0.01] * 3,
            'reason': ['-'] * 3
        }
    )
    settings = {
        'maxiter': 0, 'min_stars': 3, 'force': False, 'min_airmass_span': 0.2,
        'fallback': True, 'min_color_span': 0.3, 'rej': 3.
    }
    resid = loo(meas, 'V', ['zp', 'airmass'], {}, 0.17, settings, np.ones(3, dtype=bool))
    assert np.all(np.ma.getmaskarray(resid)), (
        'The leave-one-out residuals should be masked when the reduced fit has no degrees of '
        'freedom'
    )


def test_aperture(tmp_path, caplog):
    table, rdxdir, _ = write_night(tmp_path, wide_design())
    with caplog.at_level(logging.WARNING):
        calibrate.calibrate_all(table, rdxdir=rdxdir, mag='aper')
    assert 'without an aperture correction' in caplog.text, (
        'Calibrating aperture magnitudes should warn about the aperture correction'
    )
    content, _ = read_calibration(calibration_file(rdxdir))
    assert content['calibration']['magnitude'] == 'aper', 'The magnitude should be recorded'
    assert 'without an aperture correction' in calibration_file(rdxdir).read_text(), (
        'The calibration file should warn about the aperture correction'
    )


def test_non_ubvri(tmp_path, caplog):
    design = wide_design() + [("i'", 1.2, 4.5)]
    table, rdxdir, _ = write_night(tmp_path, design)
    with caplog.at_level(logging.WARNING):
        calibrate.calibrate_all(table, rdxdir=rdxdir, apply='image')
    assert "cannot be calibrated: calibration currently requires UBVRI" in caplog.text, (
        'A non-UBVRI filter should be skipped with a warning'
    )
    summary = Table.read(rdxdir / 'calibrate' / 'summary.ecsv')
    assert list(summary['band']) == ['B', 'V', 'R'], 'Only the UBVRI bands should be calibrated'


def test_parse_fix():
    fix = calibrate.parse_fix(['airmass=nominal', 'color=0.03', 'color2=B:0.01,V:0.02'])
    assert fix == {'airmass': 'nominal', 'color': 0.03, 'color2': {'B': 0.01, 'V': 0.02}}, (
        'The fixed terms should be parsed'
    )
    for bad in [['foo=1'], ['color'], ['color=abc']]:
        with pytest.raises(PipnickParameterError):
            calibrate.parse_fix(bad)
    with pytest.raises(PipnickParameterError):
        nightfit._resolve_fixed({'color': 'nominal'}, 'B', ['zp', 'color'], {'B': 0.3})


def test_toml_round_trip(tmp_path):
    content = {
        'calibration': {
            'setup': 's', 'night': '2026-08-20', 'magnitude': 'psf', 'input_hash': 'abc'
        },
        'band': {
            'V': {
                'status': 'ok', 'filter': 'V', 'color': 'B-V', 'terms': ['zp', 'color'],
                'coefficients': [-24.123456789, 0.1], 'errors': [0.001, 0.],
                'source': ['fit', 'fixed'], 'covariance': [[1e-6, 0.], [0., 0.]]
            }
        }
    }
    path = tmp_path / 'cal.toml'
    Calibration(content).write(path)
    parsed, edited = read_calibration(path)
    assert not edited, 'A written file should not be reported as edited'
    assert parsed['band']['V']['coefficients'] == [-24.123456789, 0.1], (
        'The coefficients should be written exactly'
    )
    assert tomllib.loads(path.read_text())['band']['V']['source'] == ['fit', 'fixed'], (
        'The file should be valid TOML'
    )
    path.write_text(path.read_text().replace('0.1]', '0.2]'))
    assert Calibration.read(path).edited, 'An edited file should be reported as edited'
