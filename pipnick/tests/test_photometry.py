"""
Tests for the photometry, using synthetic reduced images of stars with a known
Moffat point-spread function.
"""

import matplotlib.image
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from photutils.aperture import ApertureStats, CircularAperture
from photutils.psf import MoffatPSF

from pipnick.metadata.build import build_metadata
from pipnick.metadata.table import read_metadata
from pipnick.astrometry.header import apply_wcs
from pipnick.datamodel.dataset import DataSet
from pipnick.datamodel.psf import PSFModel
from pipnick.photometry import psf_photometry
from pipnick.photometry.moffat import EllipticalMoffatPSF
from pipnick.photometry.aperture_photometry import aperture_analysis
from pipnick.pipelines.photometry import photometry_all, vet_sources
from pipnick.pkg.exceptions import PipnickParameterError
from pipnick.pkg.paths import removed_table
from pipnick.pkg.version import version
from pipnick.tests.test_metadata import write_andor

# The synthetic images: 2x2-binned Andor frames (0.333 arcsec pixels), cut to
# 300x300, with a Moffat PSF of FWHM 5 pixels
SHAPE = (300, 300)
SCALE = 0.333
FWHM = 5.
BETA = 3.
ALPHA = FWHM / (2 * np.sqrt(2**(1 / BETA) - 1))
SKY = 200.
READNOISE = 10.
EXPTIME = 30.
CENTER = SkyCoord(282.77, -6.27, unit='deg')


def star_list(seed=3, nstars=25, min_dist=30.):
    """
    Return the positions and fluxes of isolated stars, away from the edges.
    """
    rng = np.random.default_rng(seed)
    x = []
    y = []
    while len(x) < nstars:
        xi, yi = rng.uniform(20, SHAPE[1] - 21, 2)
        if all(np.hypot(xi - np.array(x), yi - np.array(y)) > min_dist):
            x.append(xi)
            y.append(yi)
    return np.array(x), np.array(y), 10**rng.uniform(4.3, 5.3, nstars)


def render(x, y, flux, seed=4, error=False):
    """
    Render an image of the stars on a uniform sky, with noise.  If ``error``
    is True, also return the 1-sigma error in each pixel.
    """
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[:SHAPE[0], :SHAPE[1]]
    model = np.full(SHAPE, SKY)
    for xi, yi, fi in zip(x, y, flux):
        model += MoffatPSF(flux=fi, x_0=xi, y_0=yi, alpha=ALPHA, beta=BETA)(xx, yy)
    sigma = np.sqrt(model + READNOISE**2)
    image = model + rng.normal(0., sigma)
    return (image, sigma) if error else image


def true_wcs():
    """
    Return the true astrometric solution of the synthetic images.
    """
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.crval = [CENTER.ra.deg, CENTER.dec.deg]
    wcs.wcs.crpix = [(SHAPE[1] + 1) / 2, (SHAPE[0] + 1) / 2]
    wcs.wcs.cd = SCALE / 3600 * np.array([[-1., 0.], [0., 1.]])
    return wcs


def write_field(tmp_path, nimages=1, raw_wcs=False):
    """
    Write the metadata table and the reduced images of the synthetic stars.

    If ``raw_wcs`` is True, the images include the WCS keywords that reduced
    images inherit from the raw headers (without a pipnick solution).

    Returns the per-setup table, the reduction directory, and the stars.
    """
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    ra, dec = CENTER.to_string('hmsdms', sep=':', precision=2).split()
    for i in range(nimages):
        write_andor(
            rawdir / f'n800{i+1}.fits', OBJECT='M11 V', RA=ra, DEC=dec, XPOSURE=EXPTIME,
            AIRMASS=1.25 + 0.1 * i, **{'DATE-BEG': f'2026-08-20T04:{i:02d}:00.000'}
        )
    files, _ = build_metadata(rawdir=rawdir, filename=tmp_path / 'n8.ecsv', split=True)
    table = files[0]
    rdxdir = tmp_path / 'rdx'
    rdxdir.mkdir()
    dataset = DataSet.from_file(table, rdxdir=rdxdir)
    tag = dataset.tag
    x, y, flux = star_list()
    for product in dataset.products():
        hdr = fits.Header(
            {
                'PNSETUP': tag, 'PNRAWHSH': f'raw{product.comb_id}', 'PNTABH': 'table',
                'PROCRN': READNOISE, 'NCOMBINE': 1, 'BUNIT': 'electron'
            }
        )
        if raw_wcs:
            hdr.update(true_wcs().to_header())
        mask = np.zeros(SHAPE, dtype=np.uint8)
        mask[:,150] = 1
        image, sigma = render(x, y, flux, seed=product.comb_id, error=True)
        uncert = fits.Header({'UTYPE': 'StdDevUncertainty'})
        fits.HDUList(
            [
                fits.PrimaryHDU(data=image, header=hdr), fits.ImageHDU(data=mask, name='MASK'),
                fits.ImageHDU(data=sigma, header=uncert, name='UNCERT')
            ]
        ).writeto(product.file)
    return table, rdxdir, (x, y, flux)


def match_stars(tbl, stars, tol=2.):
    """
    Match the measured sources to the injected stars, returning the indices of
    the matched sources and stars.
    """
    x, y, _ = stars
    d = np.hypot(
        np.asarray(tbl['x'])[:,None] - x[None,:], np.asarray(tbl['y'])[:,None] - y[None,:]
    )
    near = np.argmin(d, axis=1)
    found = d[np.arange(len(tbl)), near] < tol
    return np.where(found)[0], near[found]


def test_aperture_analysis():
    data = np.ones((50, 50))
    mask = np.zeros(data.shape, dtype=bool)
    mask[10,10] = True
    flux, flux_err = aperture_analysis(
        data, [25., 10., 1.], [25., 10., 25.], 3., local_bkg=np.full(3, 0.5),
        error=np.ones_like(data), mask=mask
    )
    assert np.isclose(flux[0], 0.5 * np.pi * 9, rtol=1e-3), (
        'The flux should be the background-subtracted sum within the aperture'
    )
    expected = ApertureStats(
        data, CircularAperture([25., 25.], r=3.), error=np.ones_like(data)
    ).sum_err
    assert np.isclose(flux_err[0], expected), (
        'The flux error should be the photutils error of the aperture sum'
    )
    assert flux.mask[1] and flux.mask[2] and not flux.mask[0], (
        'Apertures with masked pixels, or beyond the image, should be masked'
    )
    _, flux_err = aperture_analysis(data, [25.], [25.], 3.)
    assert flux_err.mask[0], 'Without pixel errors, the flux errors should be masked'


def test_consolidate_groups():
    phot = Table(
        {
            'id': [1, 2, 3], 'group_id': [1, 2, 2], 'group_size': [1, 2, 2],
            'x_fit': [10., 20., 22.], 'y_fit': [5., 5., 5.], 'flux_fit': [100., 300., 100.],
            'flux_err': [1., 2., 2.], 'local_bkg': [0., 1., 1.], 'flags': [0, 1, 4]
        }
    )
    new = psf_photometry.consolidate_groups(phot)
    assert len(new) == 2, 'The two sources fit together should be consolidated'
    row = new[new['group_id'] == 2][0]
    assert row['flux_fit'] == 400. and np.isclose(row['x_fit'], 20.5), (
        'The fluxes should be summed and the positions flux-weighted'
    )
    assert row['flags'] == 5, 'The flags should be combined'
    assert len(psf_photometry.consolidate_groups(phot, preserve=[2])) == 3, (
        'Preserved groups should not be consolidated'
    )


def test_psf_model(tmp_path):
    x, y, flux = star_list()
    data = render(x, y, flux) - SKY
    for stack in ['mean', 'median']:
        fit, nstack = psf_photometry.build_moffat_psf_model(
            np.ma.MaskedArray(data), thresh=15., fwhm=5., stack=stack
        )
        fwhm = PSFModel.from_fit(fit).fwhm
        assert np.isclose(fwhm, FWHM, rtol=0.03), (
            f'The PSF FWHM from the {stack} stack should match the injected PSF'
        )
        assert nstack > 10, 'Most of the isolated stars should be used to model the PSF'
    with pytest.raises(PipnickParameterError, match='Unknown PSF fit type'):
        PSFModel('gauss', fit.par)


def test_photometry(tmp_path, caplog):
    table, rdxdir, stars = write_field(tmp_path)
    outputs = photometry_all(table, rdxdir=rdxdir, qa=True)
    assert len(outputs) == 1, 'A source table should be written for the image'
    tbl = Table.read(outputs[0])
    assert tbl.meta['pipnick_version'] == version, 'The table should record the pipnick version'
    isrc, istar = match_stars(tbl, stars)
    assert len(istar) >= 24, 'Nearly all injected stars should be found'
    assert len(tbl) - len(isrc) <= 2, 'There should be few spurious sources'

    x, y, flux = stars
    dpos = np.hypot(tbl['x'][isrc] - x[istar], tbl['y'][isrc] - y[istar])
    assert np.median(dpos) < 0.05, 'The positions should be recovered'
    ratio = np.asarray(tbl['flux_psf'][isrc]) / flux[istar]
    assert abs(np.median(ratio) - 1) < 0.01, 'The PSF fluxes should be recovered'
    pull = (np.asarray(tbl['flux_psf'][isrc]) - flux[istar]) / np.asarray(
        tbl['flux_psf_err'][isrc]
    )
    assert np.std(pull) < 2., 'The flux errors should be realistic'
    assert np.isclose(tbl.meta['psf']['fwhm_pix'], FWHM, rtol=0.03), (
        'The PSF FWHM should be recovered'
    )

    # The aperture fluxes are consistent with the PSF fluxes (scaled by the
    # fraction of the PSF within the aperture), except for masked apertures
    aper = tbl['flux_aper'][isrc]
    masked = np.absolute(np.asarray(tbl['x'][isrc]) - 150) < 2 * FWHM + 1
    assert np.all(np.ma.getmaskarray(aper)[masked]), (
        'Apertures with masked pixels should be masked'
    )
    frac = np.asarray(aper)[~masked] / flux[istar][~masked]
    assert np.std(frac) / np.mean(frac) < 0.02 and 0.85 < np.mean(frac) < 1., (
        'The aperture fluxes should be a constant fraction of the total flux'
    )

    good = np.logical_not(np.ma.getmaskarray(tbl['mag_inst']))
    assert np.allclose(
        tbl['mag_inst'][good], -2.5 * np.log10(tbl['flux_psf'][good] / EXPTIME)
    ), 'The instrumental magnitudes should be per second'
    assert np.all(np.ma.getmaskarray(tbl['ra'])) and np.all(np.ma.getmaskarray(tbl['dec'])), (
        'Without an astrometric solution, the celestial coordinates should be masked'
    )
    assert tbl.meta['airmass'] == read_metadata(table)['airmass'][0], (
        'The airmass should come from the metadata table'
    )
    assert tbl.meta['exptime'] == EXPTIME and tbl.meta['filter'] == 'V', (
        'The image metadata should be recorded'
    )
    for name in tbl.colnames:
        if tbl[name].dtype.kind == 'f':
            assert np.all(np.isfinite(np.asarray(tbl[name]))), f'{name} should not contain NaN'
    qa = rdxdir / 'photometry' / 'qa'
    assert len(sorted(qa.glob('*.png'))) == 2, 'The QA figures should be written'
    summary = Table.read(rdxdir / 'photometry' / 'summary.ecsv')
    assert summary['status'][0] == 'measured' and summary['wcs'][0] == 'none', (
        'The summary should record the measurement'
    )

    # The removed sources are listed with the reason
    removed = Table.read(removed_table(rdxdir / tbl.meta['image']))
    assert len(removed) == sum(tbl.meta['removed'].values()), (
        'Every removed source should be listed'
    )
    assert removed.meta['removed'] == tbl.meta['removed'], (
        'The removed table should record the counts of the source table'
    )
    assert set(removed['reason']) <= {'failed', 'outside', 'moved', 'nonpositive'}, (
        'Each removed source should have a reason'
    )
    for name in removed.colnames:
        if removed[name].dtype.kind == 'f':
            assert np.all(np.isfinite(np.asarray(removed[name]))), (
                f'{name} should not contain NaN'
            )



def test_vet_sources():
    phot = Table()
    phot['x_init'] = [10., 20., 30., 40., 50.]
    phot['y_init'] = [10., 20., 30., 40., 50.]
    phot['x_fit'] = [10.1, np.nan, -5., 45., 50.2]
    phot['y_fit'] = [10.1, np.nan, 30., 40., 50.1]
    phot['flux_fit'] = [100., np.nan, 50., 80., -3.]
    phot['iter_detected'] = [1, 1, 1, 2, 1]
    kept, removed, rejected = vet_sources(phot, (100, 100), 2.)
    assert list(kept['x_init']) == [10.], 'Only the credible source should be kept'
    assert removed == {'failed': 1, 'outside': 1, 'moved': 1, 'nonpositive': 1}, (
        'The removed sources should be counted by reason'
    )
    assert list(rejected['reason']) == ['failed', 'outside', 'moved', 'nonpositive'], (
        'Each removed source should be listed with its reason'
    )
    assert list(rejected['x_init']) == [20., 30., 40., 50.], (
        'The removed sources should keep their detected positions'
    )
    assert bool(rejected['x_fit'].mask[0]) and not np.any(rejected['x_fit'].mask[1:]), (
        'The failed fits should be masked, not NaN'
    )
    assert np.all(np.isfinite(rejected['flux_fit'].data)), 'The removed table should not hold NaN'


def test_photometry_qa(tmp_path):
    table, rdxdir, stars = write_field(tmp_path)
    photometry_all(table, rdxdir=rdxdir)
    qa = rdxdir / 'photometry' / 'qa'
    figures = sorted(qa.glob('*.png'))
    assert len(figures) == 2, 'By default, the QA figures should be written'
    assert all(matplotlib.image.imread(f).ndim == 3 for f in figures), (
        'The QA figures should be readable PNG images'
    )
    psf, sources = sorted(figures, key=lambda f: 'sources' in f.name)
    before = sources.stat().st_mtime_ns

    # An up-to-date image only gets the figures it lacks
    psf.unlink()
    photometry_all(table, rdxdir=rdxdir)
    assert psf.is_file(), 'A missing QA figure of an up-to-date image should be drawn'
    assert sources.stat().st_mtime_ns != before, (
        'The figures of an image are drawn together'
    )
    before = {f.name: f.stat().st_mtime_ns for f in figures}
    photometry_all(table, rdxdir=rdxdir)
    assert all(f.stat().st_mtime_ns == before[f.name] for f in figures), (
        'The figures of an up-to-date image should not be redrawn'
    )
    psf.unlink()
    photometry_all(table, rdxdir=rdxdir, qa=False)
    assert not psf.is_file(), 'No QA figure should be drawn if qa is False'


def test_coordinates_and_reruns(tmp_path, caplog):
    table, rdxdir, stars = write_field(tmp_path, raw_wcs=True)
    image = next(iter(sorted(rdxdir.glob('*.fits'))))
    ofile = photometry_all(table, rdxdir=rdxdir, qa=False)[0]
    tbl = Table.read(ofile)
    assert np.all(np.ma.getmaskarray(tbl['ra'])), (
        'The WCS inherited from the raw headers should not be used for the coordinates'
    )

    # A rerun does nothing
    before = ofile.stat().st_mtime_ns
    caplog.clear()
    photometry_all(table, rdxdir=rdxdir, qa=False)
    assert ofile.stat().st_mtime_ns == before and 'is up to date' in caplog.text, (
        'An up-to-date table should not be rewritten'
    )

    # A new astrometric solution only updates the coordinates
    wcs = true_wcs()
    apply_wcs(image, wcs)
    caplog.clear()
    photometry_all(table, rdxdir=rdxdir, qa=False)
    assert 'updated the celestial coordinates' in caplog.text, (
        'A new solution should only update the coordinates'
    )
    updated = Table.read(ofile)
    assert np.array_equal(updated['flux_psf'], tbl['flux_psf']), (
        'Updating the coordinates should not change the photometry'
    )
    expected = wcs.pixel_to_world(updated['x'], updated['y'])
    coords = SkyCoord(updated['ra'], updated['dec'], unit='deg')
    assert np.all(coords.separation(expected).arcsec < 1e-3), (
        'The coordinates should come from the solution'
    )
    assert updated.meta['wcs_method'] == 'points', 'The solution method should be recorded'

    # Changed settings, or a changed image, are measured again
    caplog.clear()
    photometry_all(table, rdxdir=rdxdir, thresh=10., qa=False)
    assert Table.read(ofile).meta['settings']['thresh'] == 10., (
        'Changed settings should trigger a new measurement'
    )
    with fits.open(image, mode='update') as hdu:
        hdu[0].header['PNTABH'] = 'edited'
    photometry_all(table, rdxdir=rdxdir, thresh=10., qa=False)
    assert Table.read(ofile).meta['image_hash'].endswith('edited'), (
        'A changed image should be measured again'
    )
    assert not np.any(np.ma.getmaskarray(Table.read(ofile)['ra'])), (
        'A new measurement should include the coordinates from the solution'
    )


def test_select_images(tmp_path):
    table, rdxdir, stars = write_field(tmp_path, nimages=2)
    outputs = photometry_all(table, rdxdir=rdxdir, images=['2'], qa=False)
    assert len(outputs) == 1 and '_C002_' in outputs[0].name, (
        'Only the selected image should be measured'
    )
    with pytest.raises(ValueError, match='No reduced science image matches'):
        photometry_all(table, rdxdir=rdxdir, images=['n9999'])


def test_elliptical_psf(tmp_path):
    table, rdxdir, stars = write_field(tmp_path)
    tbl = Table.read(photometry_all(table, rdxdir=rdxdir, fittype='ellip', qa=False)[0])
    isrc, istar = match_stars(tbl, stars)
    assert len(istar) >= 24, 'The elliptical PSF model should also find the stars'
    ratio = np.asarray(tbl['flux_psf'][isrc]) / stars[2][istar]
    assert abs(np.median(ratio) - 1) < 0.02, 'The elliptical PSF fluxes should be recovered'


def test_elliptical_moffat_psf():
    psf = EllipticalMoffatPSF(flux=1000., gamma1=3., gamma2=2., phi=0.4, alpha=3.)
    yy, xx = np.mgrid[-300:301, -300:301]
    assert np.isclose(psf(xx, yy).sum(), 1000., rtol=1e-3), (
        'The model should be normalized by its total flux'
    )
    assert [k for k, v in psf.fixed.items() if not v] == ['flux', 'x_0', 'y_0'], (
        'Only the flux and position should be free'
    )
    (ylo, yhi), (xlo, xhi) = psf.bounding_box
    assert xhi - xlo > 10 * psf.fwhm, 'The bounding box should enclose the wings of the PSF'


def test_pixel_errors(tmp_path, caplog):
    table, rdxdir, stars = write_field(tmp_path)
    image = next(iter(sorted(rdxdir.glob('*.fits'))))
    nominal = Table.read(photometry_all(table, rdxdir=rdxdir, qa=False)[0])

    # The errors come from the UNCERT extension of the reduced image
    with fits.open(image, mode='update') as hdu:
        hdu['UNCERT'].data = 2 * hdu['UNCERT'].data
    doubled = Table.read(photometry_all(table, rdxdir=rdxdir, overwrite=True, qa=False)[0])
    isrc, _ = match_stars(nominal, stars)
    jsrc, _ = match_stars(doubled, stars)
    ratio = (
        np.median(np.asarray(doubled['flux_psf_err'])[jsrc])
        / np.median(np.asarray(nominal['flux_psf_err'])[isrc])
    )
    assert np.isclose(ratio, 2., rtol=0.05), 'The flux errors should scale with the pixel errors'

    # Without pixel errors, the photometry still runs.
    # NOTE: The file must be closed before it is overwritten; Windows does not
    # allow an open file to be replaced.
    with fits.open(image, memmap=False) as hdu:
        hdus = fits.HDUList([hdu[0].copy(), hdu['MASK'].copy()])
    hdus.writeto(image, overwrite=True)
    caplog.clear()
    photometry_all(table, rdxdir=rdxdir, overwrite=True, qa=False)
    assert 'has no uncertainty' in caplog.text, 'A missing uncertainty should be reported'
