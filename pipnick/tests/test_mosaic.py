"""
Tests for the mosaics, using synthetic dithered images of stars with known
fluxes, seeing, and transparency.
"""

import logging

import numpy as np
import pytest
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import MaskedColumn, Table
from astropy.time import Time
from astropy.wcs import WCS
from photutils.aperture import ApertureStats, CircularAnnulus, CircularAperture
from photutils.psf import MoffatPSF
from scipy.signal import fftconvolve

from pipnick import metadata
from pipnick.datamodel.dataset import DataSet, check_mosaic_groups
from pipnick.datamodel.psf import PSFModel
from pipnick.mosaic import coadd, psfmatch, scaling
from pipnick.photometry.aperture_photometry import aperture_analysis
from pipnick.photometry.moffat import fwhm_to_gamma
from pipnick.pipelines import mosaic
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickParameterError
from pipnick.pkg.provenance import image_hash
from pipnick.tests.test_metadata import write_andor

SHAPE = (160, 160)
SCALE = 0.333
BETA = 2.5
SKY = 50.
READNOISE = 5.
EXPTIME = 30.
CENTER = SkyCoord(282.77, -6.27, unit='deg')


def image_wcs(dx=0., dy=0.):
    """
    Return the WCS of an image offset by (dx, dy) pixels from the center.
    """
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.crval = [CENTER.ra.deg, CENTER.dec.deg]
    wcs.wcs.crpix = [(SHAPE[1] + 1) / 2 - dx, (SHAPE[0] + 1) / 2 - dy]
    wcs.wcs.cd = SCALE / 3600 * np.array([[-1., 0.], [0., 1.]])
    return wcs


def star_field(seed=5, nstars=40):
    """
    Return the sky positions and fluxes (e-/s) of the stars, spread over the
    dithered field.
    """
    rng = np.random.default_rng(seed)
    wcs = image_wcs()
    x = rng.uniform(-40, SHAPE[1] + 40, nstars)
    y = rng.uniform(10, SHAPE[0] - 10, nstars)
    coo = wcs.pixel_to_world(x, y)
    return coo.ra.deg, coo.dec.deg, 10**rng.uniform(2.5, 3.5, nstars)


def render(wcs, stars, fwhm, transparency, background, seed):
    """
    Render an image of the stars, in electrons, with noise.  Returns the
    image, its error, and the pixel positions and fluxes of the stars in it.
    """
    rng = np.random.default_rng(seed)
    ra, dec, flux = stars
    x, y = wcs.world_to_pixel(SkyCoord(ra, dec, unit='deg'))
    yy, xx = np.mgrid[:SHAPE[0], :SHAPE[1]]
    model = np.full(SHAPE, background * EXPTIME)
    gamma = fwhm_to_gamma(fwhm, BETA)
    for xi, yi, fi in zip(x, y, flux):
        star = MoffatPSF(
            flux=fi * transparency * EXPTIME, x_0=xi, y_0=yi, alpha=gamma, beta=BETA
        )
        model += star(xx, yy)
    sigma = np.sqrt(model + READNOISE**2)
    inside = (x > 5) & (x < SHAPE[1] - 6) & (y > 5) & (y < SHAPE[0] - 6)
    return model + rng.normal(0., sigma), sigma, x[inside], y[inside], np.where(inside)[0]


def write_dithers(tmp_path, design, name='n7', stars=None):
    """
    Write the metadata table, reduced images, and photometry tables of
    dithered observations.

    Parameters
    ----------
    design : list
        For each image: the filter, the dither offset (dx, dy) in pixels, the
        FWHM (pixels), the transparency, the background (e-/s), and the UT
        hour.

    Returns
    -------
    table : :class:`pathlib.Path`
        The per-setup table, with the mosaics assigned.
    rdxdir : :class:`pathlib.Path`
        The reduction directory.
    stars : tuple
        The sky positions and fluxes of the stars.
    """
    _stars = star_field() if stars is None else stars
    rawdir = tmp_path / f'{name}_raw'
    rawdir.mkdir()
    for i, (filt, (dx, dy), _, _, _, hour) in enumerate(design):
        coo = image_wcs(dx, dy).pixel_to_world((SHAPE[1] - 1) / 2, (SHAPE[0] - 1) / 2)
        ra, dec = coo.to_string('hmsdms', sep=':', precision=2).split()
        date = Time('2026-08-20T00:00:00', scale='utc') + hour * units.hour
        write_andor(
            rawdir / f'{name}{i:03d}.fits', OBJECT=f'M11 {filt}', FILTER=filt, RA=ra, DEC=dec,
            XPOSURE=EXPTIME, AIRMASS=1.2 + 0.01 * i, **{'DATE-BEG': date.isot}
        )
    files, _ = metadata.build_metadata(
        rawdir=rawdir, filename=tmp_path / f'{name}.ecsv', split=True, assign=True,
        mosaic_radius=1.
    )
    table = files[0]
    rdxdir = tmp_path / f'{name}_rdx'
    (rdxdir / 'photometry').mkdir(parents=True)
    dataset = DataSet.from_file(table, rdxdir=rdxdir)
    tbl, tag = dataset.table, dataset.tag
    for product in dataset.products():
        i = int(str(tbl['file'][product.indx[0]])[2:5])
        filt, (dx, dy), fwhm, transparency, background, _ = design[i]
        wcs = image_wcs(dx, dy)
        image, sigma, x, y, index = render(wcs, _stars, fwhm, transparency, background, seed=i)
        hdr = fits.Header(
            {
                'PNSETUP': tag, 'PNRAWHSH': f'raw{i}', 'PNTABH': 'table', 'PNASTMTH': 'gaia',
                'PNASTWH': f'wcs{i}', 'BUNIT': 'electron'
            }
        )
        hdr.update(wcs.to_header())
        uncert = fits.Header({'UTYPE': 'StdDevUncertainty'})
        fits.HDUList(
            [
                fits.PrimaryHDU(data=image, header=hdr),
                fits.ImageHDU(data=np.zeros(SHAPE, dtype=np.uint8), name='MASK'),
                fits.ImageHDU(data=sigma, header=uncert, name='UNCERT')
            ]
        ).writeto(product.file)
        coo = wcs.pixel_to_world(x, y)
        flux = _stars[2][index] * transparency * EXPTIME
        n = len(x)
        nomask = np.zeros(n, dtype=bool)
        phot = Table()
        phot['id'] = np.arange(1, n + 1)
        phot['group_size'] = np.ones(n, dtype=int)
        phot['x'] = x
        phot['y'] = y
        phot['ra'] = MaskedColumn(coo.ra.deg, mask=nomask)
        phot['dec'] = MaskedColumn(coo.dec.deg, mask=nomask)
        phot['flux_psf'] = MaskedColumn(flux, mask=nomask)
        phot['flux_psf_err'] = MaskedColumn(np.sqrt(flux), mask=nomask)
        phot['flux_aper'] = MaskedColumn(flux, mask=nomask)
        phot['flux_aper_err'] = MaskedColumn(np.sqrt(flux), mask=nomask)
        phot['flags'] = np.zeros(n, dtype=int)
        phot.meta = {
            'image': product.file.name, 'filter': filt, 'exptime': EXPTIME,
            'psf': {
                'fittype': 'circ', 'par': list(PSFModel.circular(fwhm, BETA).par),
                'fwhm_pix': fwhm, 'fwhm_arcsec': fwhm * SCALE
            },
            'settings_hash': 'settings', 'image_hash': image_hash(hdr), 'wcs_method': 'gaia',
            'wcs_hash': f'wcs{i}'
        }
        pfile = rdxdir / 'photometry' / f'{product.file.stem}_phot.ecsv'
        phot.write(pfile, format='ascii.ecsv')
    return table, rdxdir, _stars


def three_dithers(filt='V', hour=4.):
    """
    Return the design of three dithered images with different seeing,
    transparency, and backgrounds.
    """
    return [
        (filt, (0., 0.), 4.0, 1.0, 50., hour),
        (filt, (30., 0.), 5.0, 0.8, 60., hour + 0.1),
        (filt, (60., 0.), 4.5, 0.9, 55., hour + 0.2)
    ]


def test_moffat_kernel():
    target_fwhm, target_beta = 6.0, 2.5
    size = psfmatch.kernel_size(target_fwhm, factor=10.)
    for fwhm in [4.0, 5.5]:
        model = PSFModel.circular(fwhm, 3.)
        psf = model.render(size)
        target = PSFModel.circular(target_fwhm, target_beta).render(size)
        kernel, used, ee = psfmatch.image_kernel(model, target_fwhm, target_beta)
        assert used == 'moffat', 'The Moffat kernel should be accepted'
        assert abs(ee - 1) < 0.01, 'The Moffat kernel should match the encircled energy to 1%'
        gauss = psfmatch.gaussian_kernel(fwhm, target_fwhm, size)
        ee_gauss = psfmatch.encircled_energy_ratio(psf, gauss, target, 1.5 * target_fwhm)
        assert abs(ee_gauss - 1) > abs(ee - 1), (
            'The Gaussian kernel should match the encircled energy worse than the Moffat kernel'
        )
    assert np.isclose(np.sum(kernel), 1.), 'The kernel should be normalized'
    delta = psfmatch.gaussian_kernel(6., 5., 11)
    assert delta[5,5] == 1. and np.sum(delta) == 1., (
        'A target narrower than the PSF should give a delta-function kernel'
    )


def test_normalized_convolve():
    data = np.full((60, 60), 7.)
    var = np.full((60, 60), 2.)
    good = np.ones((60, 60), dtype=bool)
    good[30,:] = False
    good[10:14,40:44] = False
    kernel = psfmatch.gaussian_kernel(2., 5., 21)
    d, v, g = psfmatch.normalized_convolve(np.where(good, data, 1e6), var, good, kernel)
    assert np.allclose(d[g], 7.), (
        'Normalized convolution should keep a constant image constant, ignoring masked pixels'
    )
    assert not np.any(g[30,:]) and not np.any(g[10:14,40:44]), (
        'Masked pixels should stay masked'
    )
    assert np.all(np.isfinite(d)) and np.all(np.isfinite(v)), 'The outputs should be finite'
    assert np.all(d[np.logical_not(g)] == 0.), 'Masked output pixels should be 0'


def test_variance_density():
    # The variance of aperture sums of convolved noise agrees with its exact
    # expectation, and the summed variance density is an upper bound on it,
    # which it approaches for apertures much larger than the kernel
    rng = np.random.default_rng(1)
    kernel = psfmatch.gaussian_kernel(1., 3., 15)
    shape = (120, 120)
    good = np.ones(shape, dtype=bool)
    var = np.full(shape, 4.)
    sums = {4.: [], 20.: []}
    for _ in range(600):
        d, v, _ = psfmatch.normalized_convolve(rng.normal(0., 2., shape), var, good, kernel)
        for radius in sums:
            sums[radius].append(float(aperture_analysis(d, [60.], [60.], radius)[0][0]))
    for radius, values in sums.items():
        aperture = CircularAperture([(60., 60.)], r=radius).to_mask(method='exact')[0]
        weights = aperture.to_image(shape)
        exact = 4. * np.sum(fftconvolve(weights, kernel[::-1,::-1], mode='same')**2)
        density = float(np.sum(weights * v))
        assert abs(np.var(values) / exact - 1) < 0.15, (
            'The variance of the aperture sums should match its exact expectation'
        )
        assert density >= exact, 'The variance density should be an upper bound'
        if radius == 20.:
            assert exact / density > 0.9, (
                'The variance density should approach the variance for large apertures'
            )


def test_solve_offsets():
    rng = np.random.default_rng(2)
    true = np.array([0., 0.1, -0.2, 0.3, 0.5])
    measurements = []
    # A chain: each image overlaps only its neighbors
    for i in range(3):
        dm = true[i] - true[i + 1] + rng.normal(0., 0.01, 10)
        measurements.append((i, i + 1, dm, np.full(10, 0.01)))
    # Image 4 overlaps with only two sources
    measurements.append((3, 4, true[3] - true[4] + np.zeros(2), np.full(2, 0.01)))
    offset, error, scaled, nused = scaling.solve_offsets(5, measurements, reference=0)
    assert np.all(np.abs(offset[:4] - true[:4]) < 0.02), (
        'The offsets should be recovered through the chain of overlapping images'
    )
    assert np.all(scaled[:4]) and not scaled[4], (
        'An image with too few common sources should not be scaled'
    )
    assert offset[4] == 0. and error[4] == 0., 'An unscaled image should have no offset'


def test_reproject_and_combine():
    stars = star_field()
    wcs_in = image_wcs()
    image, sigma, x, y, _ = render(wcs_in, stars, 4., 1., 0., seed=3)
    good = np.ones(SHAPE, dtype=bool)
    good[50:53,70:73] = False
    wcs_out, shape = coadd.output_wcs([wcs_in, image_wcs(20., 5.)], [SHAPE, SHAPE])
    bbox = coadd.footprint_bbox(wcs_in, SHAPE, wcs_out, shape)
    sub = (bbox[0].stop - bbox[0].start, bbox[1].stop - bbox[1].start)
    d, v, g = coadd.reproject_image(image, sigma**2, good, wcs_in, wcs_out[bbox], sub)
    assert np.all(np.isfinite(d)) and np.all(np.isfinite(v)), (
        'The reprojected image should be finite, with masks instead of NaN'
    )
    total_in = np.sum(image[good])
    total_out = np.sum(d[g])
    assert abs(total_out / total_in - 1) < 0.01, 'The reprojection should conserve the flux'

    # Three copies, one with a cosmic ray: rejected with three images
    images = []
    for k in range(3):
        dk = d.copy()
        if k == 1:
            dk[40,40] += 1e5
        images.append((bbox, dk, v, g))
    mean, var, coverage, used, nrej = coadd.combine(images, shape, [1., 1., 1.])
    local = (40 + bbox[0].start, 40 + bbox[1].start)
    assert nrej >= 1 and not used[1][40,40], 'The cosmic ray should be rejected'
    assert np.isclose(mean[local], d[40,40]), 'The combined value should exclude the cosmic ray'
    assert np.all(used[0][g]) and np.all(used[2][g]), (
        'The good pixels of the other images should not be rejected'
    )
    mean, _, _, used, nrej = coadd.combine(images[:2], shape, [1., 1.])
    assert nrej == 0, 'Pixels should not be rejected with fewer than three images'


def test_metadata_columns(tmp_path):
    design = three_dithers() + [('R', (0., 0.), 4., 1., 50., 5.)]
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    for i, (filt, (dx, dy), _, _, _, hour) in enumerate(design):
        coo = image_wcs(dx, dy).pixel_to_world(80, 80)
        ra, dec = coo.to_string('hmsdms', sep=':', precision=2).split()
        write_andor(
            rawdir / f'n8{i:03d}.fits', OBJECT=f'M11 {filt}', FILTER=filt, RA=ra, DEC=dec,
            **{'DATE-BEG': f'2026-08-20T0{int(hour)}:{i:02d}:00.000'}
        )
    # Another pointing, far away
    write_andor(
        rawdir / 'n8100.fits', OBJECT='SA 110 V', FILTER='V', RA='18:43:00.0', DEC='00:10:00',
        **{'DATE-BEG': '2026-08-20T06:00:00.000'}
    )
    write_andor(
        rawdir / 'n8101.fits', OBJECT='SA 110 V', FILTER='V', RA='18:43:00.5', DEC='00:10:10',
        **{'DATE-BEG': '2026-08-20T06:05:00.000'}
    )
    files, tbl = metadata.build_metadata(rawdir=rawdir, filename=tmp_path / 'n8.ecsv', split=True)
    tbl = metadata.read_metadata(files[0])
    assert np.all(tbl['mosaic'] == '') and np.all(tbl['wcs_group'] == ''), (
        'By default, no mosaics should be defined'
    )
    # NOTE: The synthetic raw frames are tiny (4x4 pixels), so the default
    # linking length of the pointing groups (half the field of view) is too
    # small, and it is given explicitly

    files, _ = metadata.build_metadata(
        rawdir=rawdir, filename=tmp_path / 'n8.ecsv', split=True, append=True, assign=True,
        mosaic_radius=1.
    )
    tbl = metadata.read_metadata(files[0])
    labels = dict(zip(tbl['file'], tbl['mosaic']))
    groups = dict(zip(tbl['file'], tbl['wcs_group']))
    assert [labels[f'n8{i:03d}.fits'] for i in range(3)] == ['M11_V'] * 3, (
        'The dithered V frames should be one mosaic, labeled without the filter in the object'
    )
    assert labels['n8003.fits'] == '', 'A single R frame should not be mosaicked'
    assert groups['n8000.fits'] == 'M11', 'The WCS group should be the pointing'
    assert labels['n8100.fits'] == 'SA_110_V' and groups['n8100.fits'] == 'SA_110', (
        'A separate pointing should be a separate WCS group'
    )

    # An edited label is kept
    edited, commented = metadata.read_metadata(files[0], include_commented=True, expand=False)
    edited['mosaic'][0] = 'mine'
    metadata.write_metadata(edited, files[0], commented=commented)
    files, _ = metadata.build_metadata(
        rawdir=rawdir, filename=tmp_path / 'n8.ecsv', split=True, append=True, assign=True,
        mosaic_radius=1.
    )
    tbl = metadata.read_metadata(files[0])
    assert tbl['mosaic'][0] == 'mine', 'An existing label should not be reassigned'


def test_check_mosaic_groups(tmp_path):
    table, rdxdir, _ = write_dithers(tmp_path, three_dithers())
    tbl = DataSet.from_file(table, rdxdir=rdxdir).table
    tbl['wcs_group'] = tbl['wcs_group'].astype('U20')
    tbl['wcs_group'][0] = 'other'
    with pytest.raises(PipnickDataModelError, match='several WCS groups'):
        check_mosaic_groups(tbl, table)


def test_mosaic(tmp_path, caplog):
    design = three_dithers('V') + three_dithers('R', hour=5.)
    # A sparse field, so that some stars are isolated
    table, rdxdir, stars = write_dithers(tmp_path, design, stars=star_field(seed=7, nstars=15))
    with caplog.at_level(logging.INFO):
        outputs = mosaic.mosaic_all(table, qa=True)
    assert sorted(f.name for f in outputs) == ['M11_R.fits', 'M11_V.fits'], (
        'There should be one mosaic per filter'
    )
    hdus = {f.stem: fits.open(f) for f in outputs}
    assert hdus['M11_V'][0].data.shape == hdus['M11_R'][0].data.shape, (
        'The mosaics of a WCS group should have the same shape'
    )
    for key in ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2']:
        assert hdus['M11_V'][0].header[key] == hdus['M11_R'][0].header[key], (
            'The mosaics of a WCS group should have the same WCS'
        )
    inputs = Table(hdus['M11_V']['INPUTS'].data)
    assert np.all(inputs['scaled']), 'All of the images should be scaled'
    # The reference is the lowest-airmass image (the first), with
    # transparency 1; the others are scaled up by 1/transparency
    expected = -2.5 * np.log10(np.array([1., 0.8, 0.9]))
    assert np.all(np.abs(inputs['offset'] - expected) < 3 * inputs['offset_err'] + 0.01), (
        'The magnitude offsets of the images should be recovered within their errors'
    )
    assert np.all(inputs['kernel'] == 'moffat'), 'The Moffat kernels should be used'
    assert hdus['M11_V'][0].header['PNTGFWHM'] == pytest.approx(1.05 * 5. * SCALE, rel=1e-3), (
        'The target FWHM should be 1.05 times the worst FWHM'
    )

    # Aperture photometry of the isolated, bright stars in the V mosaic gives
    # their fluxes on the scale of the reference image (transparency 1),
    # times the encircled energy of the target PSF
    data = hdus['M11_V'][0].data.astype(float)
    good = hdus['M11_V']['MASK'].data == 0
    assert np.all(np.isfinite(data)), 'The mosaic should be finite'
    wcs = WCS(hdus['M11_V'][0].header)
    x, y = wcs.world_to_pixel(SkyCoord(stars[0], stars[1], unit='deg'))
    target_fwhm = 1.05 * 5.
    radius = 1.5 * target_fwhm
    dist = np.hypot(x[:,None] - x[None,:], y[:,None] - y[None,:]) + np.eye(len(x)) * 1e3
    edge = 2 * target_fwhm
    use = (
        (x > edge) & (x < data.shape[1] - edge) & (y > edge) & (y < data.shape[0] - edge)
        & (np.min(dist, axis=1) > 5 * target_fwhm) & (stars[2] > 500)
    )
    annulus = CircularAnnulus(np.column_stack([x[use], y[use]]), 3 * target_fwhm, 5 * target_fwhm)
    bkg = np.asarray(ApertureStats(data, annulus, mask=np.logical_not(good)).median)
    flux, _ = aperture_analysis(
        data, x[use], y[use], radius, local_bkg=bkg, mask=np.logical_not(good)
    )
    size = psfmatch.kernel_size(target_fwhm, factor=10.)
    target = PSFModel.circular(target_fwhm, BETA).render(size)
    c = (size - 1) / 2
    yy, xx = np.mgrid[:size, :size]
    ee = np.sum(target[np.hypot(xx - c, yy - c) < radius])
    ratio = np.ma.filled(flux, 0.) / (stars[2][use] * ee)
    assert np.sum(use) >= 3, 'The test needs several isolated stars'
    assert np.median(ratio) == pytest.approx(1., abs=0.03), (
        'The mosaic should recover the star fluxes on the reference flux scale'
    )
    for h in hdus.values():
        h.close()

    # Rerun: up to date
    caplog.clear()
    with caplog.at_level(logging.INFO):
        mosaic.mosaic_all(table)
    assert caplog.text.count('is up to date') == 2, 'A rerun should not rebuild the mosaics'

    # Changing the membership of the group rebuilds all of its mosaics
    tbl, commented = metadata.read_metadata(table, include_commented=True, expand=False)
    tbl['mosaic'][np.where(tbl['mosaic'] == 'M11_R')[0][-1]] = ''
    metadata.write_metadata(tbl, table, commented=commented)
    caplog.clear()
    with caplog.at_level(logging.INFO):
        mosaic.mosaic_all(table)
    assert caplog.text.count('is up to date') == 0, (
        'A change in the membership of a WCS group should rebuild all of its mosaics'
    )


def test_mosaic_options(tmp_path, caplog):
    table, rdxdir, _ = write_dithers(tmp_path, three_dithers())
    outputs = mosaic.mosaic_all(
        table, psf_match='none', scale='none', match_background=True, outdir=tmp_path / 'out'
    )
    hdr = fits.getheader(outputs[0])
    assert hdr['PNPSFMAT'] == 'none' and 'PNTGFWHM' not in hdr, 'The PSFs should not be matched'
    inputs = Table.read(outputs[0], hdu='INPUTS')
    assert np.all(inputs['offset'] == 0.), 'The images should not be scaled'
    assert np.any(inputs['level'] != 0.), 'The background levels should be matched'
    with pytest.raises(PipnickParameterError):
        mosaic.mosaic_all(table, psf_match='bad')
    with pytest.raises(PipnickParameterError, match='output directory'):
        mosaic.mosaic_all([table, table])


def test_mixed_filters_and_nights(tmp_path, caplog):
    table, rdxdir, _ = write_dithers(tmp_path, three_dithers())
    images = mosaic.collect_images([table])
    images[1].values['filter'] = 'B'
    with caplog.at_level(logging.WARNING):
        mosaic._validate(images)
    assert 'combines images with different filters' in caplog.text, (
        'Mixed filters should give a warning'
    )
    # The images share their data set
    tbl = images[0].dataset.table
    for image in images[1:]:
        tbl['mjd'][image.product.indx] += 1.
    with pytest.raises(PipnickDataModelError, match='several nights'):
        mosaic._validate(images)


def test_two_tables(tmp_path):
    # The same mosaic observed with two setups (tables), combined into one
    stars = star_field()
    table1, _, _ = write_dithers(tmp_path, three_dithers()[:2], name='na', stars=stars)
    table2, _, _ = write_dithers(
        tmp_path, three_dithers()[2:] + three_dithers()[:1], name='nb', stars=stars
    )
    labels1 = set(metadata.read_metadata(table1)['mosaic']) - {''}
    labels2 = set(metadata.read_metadata(table2)['mosaic']) - {''}
    assert labels1 == labels2 == {'M11_V'}, 'The tables should have the same labels'
    outputs = mosaic.mosaic_all([table1, table2], outdir=tmp_path / 'joint')
    inputs = Table.read(outputs[0], hdu='INPUTS')
    assert len(inputs) == 4 and len(set(inputs['table'])) == 2, (
        'The mosaic should combine the images of both tables'
    )


def test_single_image_mosaics(tmp_path):
    # One image per filter, each its own mosaic, in one WCS group: the images
    # are only reprojected onto the WCS of the group
    design = [
        ('V', (0., 0.), 4.0, 1.0, 50., 4.), ('R', (5., 3.), 5.0, 1.0, 50., 4.1),
        ('B', (-4., 2.), 4.5, 1.0, 50., 4.2)
    ]
    table, rdxdir, _ = write_dithers(tmp_path, design)
    tbl, commented = metadata.read_metadata(table, include_commented=True, expand=False)
    for name in metadata.mosaic_columns:
        tbl[name] = tbl[name].astype('U20')
    tbl['mosaic'] = [f'M11_{f}_single' for f in tbl['filter']]
    tbl['wcs_group'] = 'M11_single'
    metadata.write_metadata(tbl, table, commented=commented)
    outputs = mosaic.mosaic_all(table)
    assert len(outputs) == 3, 'Each single-image mosaic should be constructed'
    headers = [fits.getheader(f) for f in outputs]
    for key in ['NAXIS1', 'NAXIS2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2']:
        assert len({h[key] for h in headers}) == 1, (
            'The single-image mosaics of a WCS group should have the same WCS and shape'
        )
    assert all(h['PNPSFMAT'] == 'none' and h['NCOMBINE'] == 1 for h in headers), (
        'A single image should not be convolved'
    )
    for f in outputs:
        inputs = Table.read(f, hdu='INPUTS')
        assert inputs['offset'][0] == 0. and inputs['kernel'][0] == 'none', (
            'A single image should not be scaled or convolved'
        )
        data = fits.getdata(f)
        assert np.all(np.isfinite(data)), 'The mosaic should be finite'
