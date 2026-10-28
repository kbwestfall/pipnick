"""
Tests for the astrometric calibration, using synthetic images of a synthetic
star field.  The remote services are replaced by fakes; no test contacts
astrometry.net or the Gaia archive (see ``conftest.py``).
"""

import gzip
import io
import logging
from pathlib import Path

import matplotlib.image
import numpy as np
import pytest
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from pipnick.metadata.build import build_metadata
from pipnick.metadata.table import read_metadata
from pipnick.astrometry import cache, fit, header, match, services, solve
from pipnick.astrometry.header import apply_wcs
from pipnick.astrometry.services import AstrometryNetSolver, GaiaCatalog
from pipnick.datamodel.dataset import DataSet
from pipnick.pipelines.astrometry import astrometry_all
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickError
from pipnick.pkg.version import version
from pipnick.tests.test_metadata import write_andor

# The synthetic star field
FIELD_CENTER = SkyCoord(282.77, -6.27, unit='deg')
FIELD_RADIUS = 5. / 60
NSTARS = 1500
# Images are 2x2-binned Andor frames: 0.333 arcsec pixels, cut to 400x400
SCALE = 0.333
SHAPE = (400, 400)
FWHM = 5.5
# Telescope pointing error (arcsec), the same for all images
POINTING_ERROR = (3., -2.)


def star_field(seed=1):
    """
    Return the coordinates and fluxes of the synthetic star field.
    """
    rng = np.random.default_rng(seed)
    r = FIELD_RADIUS * np.sqrt(rng.uniform(size=NSTARS))
    pa = rng.uniform(0, 360, size=NSTARS)
    coords = FIELD_CENTER.directional_offset_by(pa * units.deg, r * units.deg)
    flux = 10**rng.uniform(3., 5., size=NSTARS)
    return coords, flux


def true_wcs(center, rotation=0.):
    """
    Return the true solution of an image centered at ``center``, rotated by
    ``rotation`` degrees.
    """
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.crval = [center.ra.deg, center.dec.deg]
    wcs.wcs.crpix = [(SHAPE[1] + 1) / 2, (SHAPE[0] + 1) / 2]
    th = np.radians(rotation)
    wcs.wcs.cd = SCALE / 3600 * np.array([[-np.cos(th), np.sin(th)], [np.sin(th), np.cos(th)]])
    return wcs


def render(wcs, coords, flux, seed=2):
    """
    Render an image of the star field.
    """
    rng = np.random.default_rng(seed)
    x, y = wcs.world_to_pixel(coords)
    inside = (x > -10) & (x < SHAPE[1] + 10) & (y > -10) & (y < SHAPE[0] + 10)
    yy, xx = np.mgrid[:SHAPE[0], :SHAPE[1]]
    sig = FWHM / 2.355
    data = rng.normal(100., 5., size=SHAPE)
    for xi, yi, fi in zip(x[inside], y[inside], flux[inside]):
        near = (np.absolute(xx - xi) < 5 * sig) & (np.absolute(yy - yi) < 5 * sig)
        data[near] += fi / (2 * np.pi * sig**2) * np.exp(
            -((xx[near] - xi)**2 + (yy[near] - yi)**2) / (2 * sig**2)
        )
    return data


class FakeSolver:
    """
    Stand-in for :class:`~pipnick.astrometry.services.AstrometryNetSolver`
    that returns the true solution of the image whose stars best match the
    submitted sources.
    """
    def __init__(self, truth, coords, status=None, monitor_status=None):
        self.truth = truth
        self.coords = coords
        self.status = [] if status is None else list(status)
        self.monitor_status = [] if monitor_status is None else list(monitor_status)
        self.calls = []
        # The image for each submission ID
        self.submitted = {}
        # The headers returned, keyed by submission ID
        self.returned = {}

    def _identify(self, x, y):
        best = None
        best_dist = np.inf
        for name, wcs in self.truth.items():
            sx, sy = wcs.world_to_pixel(self.coords)
            d = np.min(np.hypot(sx[None,:] - x[:,None], sy[None,:] - y[:,None]), axis=1)
            if np.median(d) < best_dist:
                best_dist = np.median(d)
                best = name
        return best

    def header(self, name):
        """
        Return the true solution as astrometry.net would: with ``EQUINOX =
        2000.0`` and no ``RADESYS``, i.e., labeled as FK5.
        """
        hdr = fits.Header({'SIMPLE': True, 'BITPIX': 8, 'NAXIS': 0})
        hdr.update(self.truth[name].to_header())
        hdr.remove('RADESYS', ignore_missing=True)
        hdr['EQUINOX'] = 2000.
        hdr['COMMENT'] = 'Synthetic astrometry.net solution'
        return hdr

    def _result(self, name):
        submission_id = f'sub{len(self.calls)}'
        self.submitted[submission_id] = name
        status = self.status.pop(0) if len(self.status) > 0 else 'solved'
        if status != 'solved':
            return status, None, submission_id
        self.returned[submission_id] = self.header(name)
        return 'solved', self.returned[submission_id].copy(), submission_id

    def solve_sources(self, x, y, width, height, hints):
        name = self._identify(np.asarray(x), np.asarray(y))
        self.calls.append(('sources', name, hints))
        return self._result(name)

    def solve_image(self, path, hints):
        # NOTE: Memory mapping would keep the file open, and Windows does not
        # allow the pipeline to delete it
        data = fits.getdata(path, memmap=False)
        self.calls.append(('image', path, data))
        return self._result(self.uploaded)

    def monitor(self, submission_id):
        self.calls.append(('monitor', submission_id))
        status = self.monitor_status.pop(0) if len(self.monitor_status) > 0 else 'solved'
        if status != 'solved':
            return status, None, submission_id
        return 'solved', self.header(self.submitted[submission_id]), submission_id


class FakeCatalog:
    """
    Stand-in for :class:`~pipnick.astrometry.services.GaiaCatalog` that
    returns the synthetic star field.
    """
    def __init__(self, coords, flux):
        self.coords = coords
        self.mag = -2.5 * np.log10(flux)
        self.calls = 0

    def sources(self, center, radius, epoch):
        self.calls += 1
        inside = self.coords.separation(center).deg <= radius
        return self.coords[inside], self.mag[inside]


# Image centers as offsets from the field center (arcsec) and rotations (deg):
# three images within 1 arcmin of each other, and one ~3 arcmin away.
IMAGES = {
    'n7001': ((0., 0.), 0.), 'n7002': ((20., -15.), 0.), 'n7003': ((-40., 30.), 0.05),
    'n7004': ((150., 90.), 0.)
}


def write_field(tmp_path, images=IMAGES, stale=None):
    """
    Write the metadata table and the reduced images of the synthetic field.

    ``stale`` maps image names to the name of the image whose pointing they
    are given, to mimic a stale pointing in the raw header.

    Returns the per-setup table, the reduction directory, the true solutions
    keyed by the reduced image name, and the star field.
    """
    coords, flux = star_field()
    rawdir = tmp_path / 'raw'
    rawdir.mkdir()
    centers = {}
    pointings = {}
    for name, (offset, rotation) in images.items():
        center = FIELD_CENTER.spherical_offsets_by(
            offset[0] * units.arcsec, offset[1] * units.arcsec
        )
        pointings[name] = center.spherical_offsets_by(
            POINTING_ERROR[0] * units.arcsec, POINTING_ERROR[1] * units.arcsec
        )
        centers[f'{name}.fits'] = (center, rotation)
    for i, name in enumerate(images.keys()):
        pointing = pointings[name if stale is None else stale.get(name, name)]
        ra, dec = pointing.to_string('hmsdms', sep=':', precision=2).split()
        write_andor(
            rawdir / f'{name}.fits', OBJECT='M11 V', RA=ra, DEC=dec,
            **{'DATE-BEG': f'2026-08-20T04:{i:02d}:00.000'}
        )
    files, _ = build_metadata(rawdir=rawdir, filename=tmp_path / 'n7.ecsv', split=True)
    table = files[0]
    rdxdir = tmp_path / 'rdx'
    rdxdir.mkdir()
    dataset = DataSet.from_file(table, rdxdir=rdxdir)
    tbl, tag = dataset.table, dataset.tag
    truth = {}
    for product in dataset.products():
        center, rotation = centers[tbl['file'][product.indx[0]]]
        wcs = true_wcs(center, rotation)
        hdr = fits.Header({'PNSETUP': tag, 'OBJECT': 'M11'})
        mask = np.zeros(SHAPE, dtype=np.uint8)
        mask[:,100] = 1
        fits.HDUList(
            [
                fits.PrimaryHDU(data=render(wcs, coords, flux), header=hdr),
                fits.ImageHDU(data=mask, name='MASK')
            ]
        ).writeto(product.file)
        truth[product.file.name] = wcs
    return table, rdxdir, truth, coords, flux


def solution_error(path, wcs):
    """
    Return the maximum difference in pixels between the solution in an image
    and its true solution, over a grid of positions.
    """
    solved = fit.wcs_from_header(fits.getheader(path))
    gx, gy = np.meshgrid(np.linspace(0, SHAPE[1] - 1, 5), np.linspace(0, SHAPE[0] - 1, 5))
    px, py = wcs.world_to_pixel(solved.pixel_to_world(gx.ravel(), gy.ravel()))
    return float(np.max(np.hypot(px - gx.ravel(), py - gy.ravel())))


def test_fit_wcs():
    wcs = true_wcs(FIELD_CENTER, rotation=0.3)
    rng = np.random.default_rng(3)
    x = rng.uniform(0, SHAPE[1] - 1, 40)
    y = rng.uniform(0, SHAPE[0] - 1, 40)
    coords = wcs.pixel_to_world(x, y)
    solved, stats = fit.fit_wcs(x, y, coords)
    assert stats['npts'] == 40 and stats['rms'] < 1e-3, 'Exact points should be fit exactly'
    px, py = wcs.world_to_pixel(solved.pixel_to_world(x, y))
    assert np.max(np.hypot(px - x, py - y)) < 1e-3, (
        'The fit should use 0-based pixel coordinates, independent of the astropy version'
    )

    # Noise and outliers
    nx = x + rng.normal(0, 0.1, 40)
    ny = y + rng.normal(0, 0.1, 40)
    nx[:3] += 20.
    solved, stats = fit.fit_wcs(nx, ny, coords)
    assert not np.any(stats['used'][:3]) and np.all(stats['used'][3:]), (
        'Outliers, and only outliers, should be rejected'
    )
    assert stats['rms'] < 0.1, 'The RMS residual should reflect the noise'

    with pytest.raises(PipnickDataModelError, match='At least 3'):
        fit.fit_wcs(x[:2], y[:2], coords[:2])


def test_fit_wcs_origin(monkeypatch):
    # Versions of astropy with the origin argument must be told the pixels are 0-based
    received = {}

    def fake_fit(
        xy, world_coords, origin=None, proj_point='center', projection='TAN', sip_degree=None
    ):
        received['origin'] = origin
        return true_wcs(FIELD_CENTER)

    monkeypatch.setattr(fit, 'fit_wcs_from_points', fake_fit)
    fit._fit_wcs_from_points(np.zeros(3), np.zeros(3), SkyCoord([0, 1, 2], [0, 0, 0], unit='deg'))
    assert received['origin'] == 0, 'The origin should be 0 when astropy accepts it'


def test_match_frames():
    rng = np.random.default_rng(4)
    xref = rng.uniform(0, 400, 60)
    yref = rng.uniform(0, 400, 60)
    rotation = np.radians(0.2)
    shift = np.array([-35.2, 22.7])
    # The reference positions are q = R p + t, so p = R^-1 (q - t)
    c = np.cos(-rotation)
    s = np.sin(-rotation)
    x = c * (xref - shift[0]) - s * (yref - shift[1])
    y = s * (xref - shift[0]) + c * (yref - shift[1])
    # Keep a subset, and add spurious sources
    x = np.append(x[:40], rng.uniform(0, 400, 15))
    y = np.append(y[:40], rng.uniform(0, 400, 15))
    result = match.match_frames(x, y, xref, yref, expected=(-30., 20.), search=20.)
    assert np.isclose(result['rotation'], rotation, atol=1e-4), 'The rotation should be found'
    assert np.allclose(result['shift'], shift, atol=0.01), 'The shift should be found'
    common = result['indx'] < 40
    assert np.array_equal(result['indx_ref'][common], result['indx'][common]), (
        'Sources in both sets should be matched to the correct partners'
    )
    assert np.sum(common) >= 38, 'Most common sources should be matched'
    assert np.sum(np.logical_not(common)) <= 2, (
        'Spurious sources should only rarely coincide with a reference source by chance'
    )


def test_propagate_proper_motion():
    catalog = Table(
        {
            'ra': [10., 10.], 'dec': [60., 60.], 'pmra': np.ma.MaskedArray([0., 0.]),
            'pmdec': np.ma.MaskedArray([1000., 0.], mask=[False, True]),
            'ref_epoch': [2016., 2016.]
        }
    )
    coords = services.propagate_proper_motion(catalog, 2026.)
    assert np.isclose(coords[0].dec.deg - 60., 10. / 3600), (
        'A star should move by its proper motion times the elapsed time'
    )
    assert coords[1].dec.deg == 60., 'A star with a masked proper motion should not move'


def test_icrs_header():
    # As returned by astrometry.net, and as written by astropy after reading it
    hdr = true_wcs(FIELD_CENTER).to_header()
    hdr.remove('RADESYS', ignore_missing=True)
    hdr['EQUINOX'] = 2000.
    assert fit.wcs_from_header(hdr).wcs.radesys == 'FK5', (
        'astropy should read an astrometry.net solution as FK5'
    )
    for h in [hdr, fit.wcs_from_header(hdr).to_header()]:
        relabeled = services.icrs_header(h)
        assert relabeled['RADESYS'] == 'ICRS' and 'EQUINOX' not in relabeled, (
            'The solution should be relabeled as ICRS'
        )
        wcs = fit.wcs_from_header(relabeled)
        assert wcs.wcs.radesys == 'ICRS', 'astropy should read the relabeled solution as ICRS'
        assert np.allclose(wcs.wcs.crval, true_wcs(FIELD_CENTER).wcs.crval), (
            'Relabeling should not change the solution'
        )
    assert 'EQUINOX' in hdr, 'The input header should not be modified'


def test_image_center_frame():
    wcs = true_wcs(FIELD_CENTER)
    wcs.wcs.radesys = 'FK5'
    wcs.wcs.equinox = 2000.
    center = fit.image_center(wcs, SHAPE)
    assert center.frame.name == 'icrs', 'Image centers should be returned in the ICRS'
    assert center.is_equivalent_frame(FIELD_CENTER), (
        'Converted centers should carry no FK5 attributes, so they can be compared with ICRS '
        'coordinates'
    )
    center.spherical_offsets_to(FIELD_CENTER)
    assert center.separation(FIELD_CENTER).arcsec < 0.1, (
        'The FK5 center should be converted, not relabeled'
    )


def test_api_key(monkeypatch, caplog):
    monkeypatch.delenv('ASTROMETRY_NET_API_KEY', raising=False)
    monkeypatch.setattr(services.astrometry_net.conf, 'api_key', '')
    with pytest.raises(PipnickError, match='No astrometry.net API key'):
        AstrometryNetSolver()
    monkeypatch.setenv('ASTROMETRY_NET_API_KEY', 'not-a-real-key')
    AstrometryNetSolver()
    assert 'not-a-real-key' not in caplog.text, 'The API key should never be logged'


def test_blocked_services(tmp_path):
    # The test configuration blocks any contact with the remote services
    solver = AstrometryNetSolver(api_key='not-a-real-key')
    with pytest.raises(RuntimeError, match='must not contact'):
        solver.solve_sources(np.zeros(3), np.zeros(3), 10, 10, {})
    with pytest.raises(RuntimeError, match='must not contact'):
        GaiaCatalog(tmp_path).sources(FIELD_CENTER, 0.1, 2026.)


def test_astrometry(tmp_path, caplog):
    table, rdxdir, truth, coords, flux = write_field(tmp_path)
    solver = FakeSolver(truth, coords)
    catalog = FakeCatalog(coords, flux)
    solved = astrometry_all(table, rdxdir=rdxdir, solver=solver, catalog=catalog, qa=False)

    assert len(solved) == 4, 'All images should be solved'
    assert len(solver.calls) == 1, 'Only one image should be submitted'
    assert solver.calls[0][2]['publicly_visible'] == 'n', 'Solutions should not be public'
    methods = {p.name: fits.getheader(p)['PNASTMTH'] for p in solved}
    submitted = solver.calls[0][1]
    assert methods[submitted] == 'astrometry.net', 'The submitted image should be marked'
    assert sorted(methods.values()) == ['astrometry.net', 'bootstrap', 'bootstrap', 'gaia'], (
        'The images near the submitted one should be bootstrapped, and the distant one solved '
        'with the catalog'
    )
    # NOTE: An image can be bootstrapped from a bootstrapped image, so the
    # errors can accumulate; 0.2 pixel is 0.07 arcsec
    for p in solved:
        assert solution_error(p, truth[p.name]) < 0.2, (
            f'The solution for {p.name} should match the truth to better than 0.2 pixel'
        )
        hdr = fits.getheader(p)
        assert hdr['RADESYS'] == 'ICRS' and 'EQUINOX' not in hdr, (
            f'The solution for {p.name} should be in the ICRS'
        )
    assert catalog.calls == 1, 'Only the distant image should need the catalog'
    saved = rdxdir / 'astrometry' / f'anet_{fits.getheader(rdxdir / submitted)["PNASTIN"]}.fits'
    assert (
        fits.Header.fromfile(saved).tostring() == next(iter(solver.returned.values())).tostring()
    ), 'The astrometry.net solution should be saved exactly as returned'
    index = Table.read(rdxdir / 'astrometry' / 'astrometry_net.ecsv')
    assert list(index['status']) == ['solved'] and index['solution'][0] == saved.name, (
        'The submission should be recorded in the index'
    )

    # Rerunning changes nothing
    before = {p.name: p.stat().st_mtime_ns for p in solved}
    astrometry_all(table, rdxdir=rdxdir, solver=solver, catalog=catalog, qa=False)
    assert len(solver.calls) == 1, 'A rerun should not submit anything'
    assert {p.name: p.stat().st_mtime_ns for p in solved} == before, (
        'A rerun should not change solved images'
    )

    # Rewriting the images (e.g., by pipnick_reduce) restores the cached
    # solution and redoes the rest, without submitting
    for p in solved:
        with fits.open(p, mode='update') as hdu:
            header._remove_wcs(hdu[0].header)
    caplog.clear()
    astrometry_all(table, rdxdir=rdxdir, solver=solver, catalog=catalog, qa=False)
    assert len(solver.calls) == 1, 'Rewritten images should not be submitted again'
    assert 'cached astrometry.net solution' in caplog.text, 'The cache should be used'
    for p in solved:
        assert solution_error(p, truth[p.name]) < 0.2, f'{p.name} should be solved again'

    # An astrometry.net solution applied with the FK5 label is relabeled from
    # the cache, without submitting
    submitted_path = rdxdir / submitted
    with fits.open(submitted_path, mode='update') as hdu:
        hdu[0].header['RADESYS'] = 'FK5'
        hdu[0].header['EQUINOX'] = 2000.
    astrometry_all(table, rdxdir=rdxdir, solver=solver, catalog=catalog, qa=False)
    hdr = fits.getheader(submitted_path)
    assert hdr['RADESYS'] == 'ICRS' and 'EQUINOX' not in hdr, (
        'An FK5-labeled astrometry.net solution should be relabeled as ICRS'
    )
    assert len(solver.calls) == 1, 'Relabeling should not submit anything'

    # Re-solving submits again
    astrometry_all(table, rdxdir=rdxdir, solver=solver, catalog=catalog, resolve=True, qa=False)
    assert len(solver.calls) == 2, 'Resolve should submit again'


def test_reference_changed(tmp_path, caplog):
    table, rdxdir, truth, coords, flux = write_field(tmp_path)
    solver = FakeSolver(truth, coords)
    astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), qa=False
    )
    # NOTE: The order of glob is platform-dependent, so the files are sorted
    bootstrapped = [
        p for p in sorted(rdxdir.glob('*.fits')) if fits.getheader(p)['PNASTMTH'] == 'bootstrap'
    ]
    ref = rdxdir / fits.getheader(bootstrapped[0])['PNASTREF']
    # The reference may itself be bootstrapped; it gets a new (user) solution
    dependents = [p for p in bootstrapped if p != ref]
    assert len(dependents) > 0, 'At least one image should depend on the changed reference'

    # A new solution for the reference makes its dependents out of date
    wcs = fit.wcs_from_header(fits.getheader(ref))
    wcs.wcs.crval = wcs.wcs.crval + np.array([0., 1e-5])
    apply_wcs(ref, wcs)
    caplog.clear()
    astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), qa=False
    )
    assert 'reference solution for' in caplog.text, 'The changed reference should be noticed'
    for p in dependents:
        hdr = fits.getheader(p)
        assert hdr['PNASTRH'] == fits.getheader(rdxdir / hdr['PNASTREF'])['PNASTWH'], (
            'The dependents should be solved again from the current reference solution'
        )
    assert len(solver.calls) == 1, 'Redoing the dependents should not submit anything'


def test_points(tmp_path):
    table, rdxdir, truth, coords, flux = write_field(tmp_path, images={'n7001': ((0., 0.), 0.)})
    image = next(rdxdir.glob('*.fits'))
    wcs = truth[image.name]
    rng = np.random.default_rng(5)
    x = rng.uniform(0, SHAPE[1] - 1, 10)
    y = rng.uniform(0, SHAPE[0] - 1, 10)
    solution, stats = fit.fit_wcs(x, y, wcs.pixel_to_world(x, y))
    apply_wcs(image, solution, stats=stats)
    hdr = fits.getheader(image)
    assert hdr['PNASTMTH'] == 'points' and hdr['PNASTNPT'] == 10, (
        'The user solution and its statistics should be recorded'
    )
    assert hdr['PNASTVER'] == version, 'The solution should record the pipnick version'
    assert (rdxdir / 'astrometry' / f'{image.stem}_points.fits').is_file(), (
        'The user solution should be saved'
    )

    # The user solution is restored after the image is rewritten, without submitting
    with fits.open(image, mode='update') as hdu:
        header._remove_wcs(hdu[0].header)
    solver = FakeSolver(truth, coords)
    astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), qa=False
    )
    assert fits.getheader(image)['PNASTMTH'] == 'points', 'The user solution should be restored'
    assert len(solver.calls) == 0, 'Nothing should be submitted'
    assert solution_error(image, wcs) < 1e-3, 'The restored solution should be unchanged'


def test_pending_and_failed(tmp_path):
    images = {'n7001': ((0., 0.), 0.), 'n7004': ((150., 90.), 0.)}
    table, rdxdir, truth, coords, flux = write_field(tmp_path, images=images)

    # The first submission is pending and finishes later; the other image is
    # solved from the catalog once the first is solved
    solver = FakeSolver(truth, coords, status=['pending'])
    solved = astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), qa=False
    )
    assert [c[0] for c in solver.calls] == ['sources', 'monitor'], (
        'The pending submission should be monitored'
    )
    assert len(solved) == 2, 'Both images should be solved'
    for p in solved:
        assert solution_error(p, truth[p.name]) < 0.1, (
            'The monitored solution should be applied to the image that was submitted'
        )

    # A submission still pending at the end of a run is monitored by the next
    (tmp_path / 'c').mkdir()
    table, rdxdir, truth, coords, flux = write_field(tmp_path / 'c', images=images)
    solver = FakeSolver(truth, coords, status=['pending'], monitor_status=['pending'])
    solved = astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), qa=False
    )
    assert len(solved) == 0 and [c[0] for c in solver.calls] == ['sources', 'monitor'], (
        'Nothing should be solved, or submitted, while a submission is pending'
    )
    solved = astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), qa=False
    )
    assert [c[0] for c in solver.calls] == ['sources', 'monitor', 'monitor'], (
        'The pending submission should be monitored again, not resubmitted'
    )
    assert len(solved) == 2, 'Both images should be solved once the submission finishes'

    # A failed submission is not repeated
    (tmp_path / 'b').mkdir()
    table, rdxdir, truth, coords, flux = write_field(tmp_path / 'b', images=images)
    solver = FakeSolver(truth, coords, status=['failed'])
    solved = astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), max_submissions=1,
        qa=False
    )
    assert len(solved) == 0, 'Nothing can be solved without a solved image'
    failed = solver.calls[0][1]
    astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), qa=False
    )
    assert [c[1] for c in solver.calls].count(failed) == 1, (
        'An image whose submission failed should not be resubmitted'
    )


def test_max_submissions(tmp_path):
    table, rdxdir, truth, coords, flux = write_field(tmp_path)
    solver = FakeSolver(truth, coords)
    solved = astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), max_submissions=0,
        qa=False
    )
    assert len(solver.calls) == 0 and len(solved) == 0, 'The submission limit should be respected'


def test_dry_run(tmp_path, caplog):
    table, rdxdir, truth, coords, flux = write_field(tmp_path)
    # The default solver and catalog would contact the services if used
    solved = astrometry_all(table, rdxdir=rdxdir, dry_run=True, qa=False)
    assert len(solved) == 0, 'A dry run should not solve anything'
    assert 'At least 1 and' in caplog.text, 'The number of submissions should be reported'
    assert 'Group 2' in caplog.text, 'The pointing groups should be reported'
    written = sorted((rdxdir / 'astrometry' / 'dry_run').glob('*_sources.ecsv'))
    assert len(written) == 4, 'The source list of each image should be written'
    assert 'center_ra' in Table.read(written[0]).meta['hints'], 'The hints should be written'
    for p in rdxdir.glob('*.fits'):
        assert 'PNASTMTH' not in fits.getheader(p), 'A dry run should not change the images'


def test_upload_image(tmp_path):
    images = {'n7001': ((0., 0.), 0.)}
    table, rdxdir, truth, coords, flux = write_field(tmp_path, images=images)
    solver = FakeSolver(truth, coords)
    solver.uploaded = next(iter(truth.keys()))
    solved = astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), upload_image=True,
        qa=False
    )
    assert len(solved) == 1 and solver.calls[0][0] == 'image', 'The image should be uploaded'
    assert np.all(solver.calls[0][2][:,100] == 0), 'Masked pixels should be set to 0'
    assert not solver.calls[0][1].exists(), 'The uploaded copy should be removed'


# Two images with the stale pointing of the first, ~3.3 arcmin from their true
# centers
STALE_IMAGES = {
    'n7001': ((0., 0.), 0.), 'n7002': ((170., 100.), 0.), 'n7003': ((175., 95.), 0.)
}
STALE = {'n7002': 'n7001', 'n7003': 'n7001'}


def solve_first(rdxdir, truth):
    """
    Apply the true solution to the first image, as provided by the user.
    """
    first = sorted(truth.keys())[0]
    apply_wcs(rdxdir / first, truth[first])
    return first


def test_wide_search(tmp_path, caplog):
    table, rdxdir, truth, coords, flux = write_field(
        tmp_path, images=STALE_IMAGES, stale=STALE
    )
    first, stale_1, stale_2 = sorted(truth.keys())
    solve_first(rdxdir, truth)
    solver = FakeSolver(truth, coords)
    solved = astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), qa=False
    )

    assert len(solved) == 3 and len(solver.calls) == 0, (
        'Images with a stale pointing should be solved without submitting'
    )
    hdr = fits.getheader(rdxdir / stale_1)
    assert hdr['PNASTMTH'] == 'gaia', 'The first stale image should be solved by the wide search'
    tbl = read_metadata(table)
    pointing = SkyCoord(tbl['ra'][1], tbl['dec'][1], unit=('hourangle', 'deg'))
    center = fit.image_center(truth[stale_1], SHAPE)
    assert np.isclose(hdr['PNASTPOE'], center.separation(pointing).arcmin, atol=0.01), (
        'The pointing error should be recorded'
    )
    assert 'consider correcting the RA and Dec in the metadata table' in caplog.text, (
        'The pointing error should be reported'
    )
    hdr = fits.getheader(rdxdir / stale_2)
    assert hdr['PNASTMTH'] == 'bootstrap' and hdr['PNASTREF'] == stale_1, (
        'An image with the same stale pointing should use the measured pointing error, and be '
        'bootstrapped from the first'
    )
    for p in solved:
        assert solution_error(p, truth[p.name]) < 0.1, f'{p.name} should be solved correctly'

    # A rerun keeps the solutions
    before = {p.name: p.stat().st_mtime_ns for p in solved}
    astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), qa=False
    )
    assert {p.name: p.stat().st_mtime_ns for p in solved} == before, (
        'A rerun should not change the solutions'
    )

    # If the second stale image is rewritten, it is solved again using the
    # pointing error of the first, without the wide search
    with fits.open(rdxdir / stale_2, mode='update') as hdu:
        header._remove_wcs(hdu[0].header)
    caplog.clear()
    catalog = FakeCatalog(coords, flux)
    astrometry_all(table, rdxdir=rdxdir, solver=solver, catalog=catalog, qa=False)
    assert catalog.calls == 0, 'The restored pointing error should be used, without the catalog'
    assert f'{stale_2} is centered at' in caplog.text, (
        'The incorrect pointing should be reported for every affected image'
    )
    assert 'PNASTPOE' in fits.getheader(rdxdir / stale_2), (
        'The pointing error should be recorded for every affected image'
    )
    assert fits.getheader(rdxdir / stale_2)['PNASTREF'] == stale_1, (
        'The image should be bootstrapped from the first stale image'
    )


def test_wide_search_disabled(tmp_path):
    table, rdxdir, truth, coords, flux = write_field(
        tmp_path, images=STALE_IMAGES, stale=STALE
    )
    solve_first(rdxdir, truth)
    solver = FakeSolver(truth, coords)
    solved = astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), wide_search=0.,
        max_submissions=0, qa=False
    )
    assert len(solved) == 1, 'Without the wide search, the stale images should not be solved'


def test_invalidate_chains():
    # C was bootstrapped from B, which was bootstrapped from A
    images = []
    chain = [
        ('a', 'astrometry.net', None), ('b', 'bootstrap', 'a.fits'), ('c', 'bootstrap', 'b.fits'),
        ('d', 'bootstrap', 'a.fits')
    ]
    for name, method, ref in chain:
        hdr = fits.Header({'PNASTMTH': method})
        if ref is not None:
            hdr['PNASTREF'] = ref
        image = solve.AstrometryImage(
            Path(f'{name}.fits'), hdr, SHAPE, FIELD_CENTER, 2026., SCALE, FWHM, None, {}, ''
        )
        image.status = 'solved'
        image.method = method
        images.append(image)
    by_name = {image.name: image for image in images}

    cache.invalidate_chains(images, by_name)
    assert all(image.status == 'solved' for image in images), (
        'Valid chains should not be invalidated'
    )
    by_name['b.fits'].status = None
    cache.invalidate_chains(images, by_name)
    assert [image.status for image in images] == ['solved', None, None, 'solved'], (
        'Invalidating an image should invalidate the chain bootstrapped from it, and only that '
        'chain'
    )


def fake_gaia_query(catalog, query, output_file):
    """
    Stand-in for :meth:`~pipnick.astrometry.services.GaiaCatalog._query` that
    writes the synthetic star field as a gzip-compressed VOTable, as the
    archive would.
    """
    coords, flux = star_field()
    catalog.nqueries = getattr(catalog, 'nqueries', 0) + 1
    buffer = io.BytesIO()
    Table(
        {
            'ra': coords.ra.deg, 'dec': coords.dec.deg, 'pmra': np.zeros(len(coords)),
            'pmdec': np.full(len(coords), 100.), 'ref_epoch': np.full(len(coords), 2016.),
            'phot_g_mean_mag': -2.5 * np.log10(flux) + 25.
        }
    ).write(buffer, format='votable')
    with gzip.open(output_file, 'wb') as f:
        f.write(buffer.getvalue())


def test_gaia_catalog(tmp_path, monkeypatch):
    monkeypatch.setattr(GaiaCatalog, '_query', fake_gaia_query)
    catalog = GaiaCatalog(tmp_path / 'astrometry')
    coords, mag = catalog.sources(FIELD_CENTER, 2. / 60, 2026.)
    saved = tmp_path / 'astrometry' / 'gaia_0000.vot.gz'
    assert saved.is_file() and catalog.nqueries == 1, 'The query result should be saved'
    content = saved.read_bytes()
    assert np.all(coords.separation(FIELD_CENTER).arcmin <= 2.), (
        'Only sources in the requested cone should be returned'
    )
    # 100 mas/yr over 10 yr moves each star 1 arcsec north; moving the returned
    # positions back should recover the catalog positions exactly
    raw = Table.read(saved, format='votable')
    back = coords.spherical_offsets_by(0. * units.arcsec, -1. * units.arcsec)
    _, sep, _ = back.match_to_catalog_sky(SkyCoord(raw['ra'], raw['dec'], unit='deg'))
    assert np.all(sep.arcsec < 1e-6), 'Positions should be propagated to the epoch'

    # A covered cone reuses the saved query, which is not modified
    catalog.sources(
        FIELD_CENTER.directional_offset_by(0. * units.deg, 0.5 * units.arcmin), 1. / 60, 2026.
    )
    assert catalog.nqueries == 1, 'A covered cone should not be queried again'
    assert saved.read_bytes() == content, 'The saved query result should never be modified'
    index = Table.read(tmp_path / 'astrometry' / 'gaia_queries.ecsv')
    assert len(index) == 1 and 'gaiadr3.gaia_source' in index['query'][0], (
        'The query should be recorded in the index'
    )

    # A deeper catalog needs a new query
    GaiaCatalog(tmp_path / 'astrometry', mag_limit=20.).sources(FIELD_CENTER, 1. / 60, 2026.)
    assert (tmp_path / 'astrometry' / 'gaia_0001.vot.gz').is_file(), (
        'A query that is not deep enough should not be reused'
    )


def test_select_images(tmp_path, caplog):
    table, rdxdir, truth, coords, flux = write_field(tmp_path)
    solver = FakeSolver(truth, coords)
    catalog = FakeCatalog(coords, flux)
    astrometry_all(table, rdxdir=rdxdir, solver=solver, catalog=catalog, qa=False)
    names = sorted(truth.keys())
    for name in names[1:3]:
        with fits.open(rdxdir / name, mode='update') as hdu:
            header._remove_wcs(hdu[0].header)

    # Only the selected image (by comb_id) is solved
    solved = astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=catalog, images=['2'], qa=False
    )
    assert names[1] in [p.name for p in solved] and names[2] not in [p.name for p in solved], (
        'Only the selected image should be solved'
    )
    summary = Table.read(rdxdir / 'astrometry' / 'summary.ecsv')
    assert summary['status'][summary['image'] == names[2]][0] == 'skipped', (
        'Images that are not selected should be reported as skipped'
    )

    # Re-solving a selected image changes only that image
    submitted = solver.calls[0][1]
    before = {name: (rdxdir / name).stat().st_mtime_ns for name in names}
    astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=catalog, images=[submitted], resolve=True,
        max_submissions=1, qa=False
    )
    changed = [name for name in names if (rdxdir / name).stat().st_mtime_ns != before[name]]
    assert changed == [submitted], 'Resolve should only change the selected image'
    assert len(solver.calls) <= 2, 'Resolve should submit at most the selected image'

    with pytest.raises(ValueError, match='No reduced science image matches'):
        astrometry_all(table, rdxdir=rdxdir, solver=solver, catalog=catalog, images=['n9999'])


def test_summary(tmp_path, caplog):
    table, rdxdir, truth, coords, flux = write_field(tmp_path)
    solver = FakeSolver(truth, coords)
    astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), max_submissions=0,
        qa=False
    )
    summary = Table.read(rdxdir / 'astrometry' / 'summary.ecsv')
    assert list(summary['status']) == ['unsolved'] * 4, (
        'Nothing can be solved without a submission'
    )
    assert all(r != '-' for r in summary['reason']), (
        'The reason each image was not solved should be recorded'
    )
    assert any('--max_submissions' in s for s in summary['suggestion']), (
        'The summary should suggest allowing a submission'
    )
    assert 'is not solved' in caplog.text, 'Unsolved images should be reported'

    astrometry_all(
        table, rdxdir=rdxdir, solver=solver, catalog=FakeCatalog(coords, flux), qa=False
    )
    summary = Table.read(rdxdir / 'astrometry' / 'summary.ecsv')
    assert list(summary['status']) == ['solved'] * 4, 'All images should be solved'
    assert np.all(summary['rms'][summary['method'] != 'astrometry.net'] < 0.3), (
        'The summary should give the RMS of the solutions'
    )
    assert summary['nmatch'].mask[summary['method'] == 'astrometry.net'][0], (
        'Unavailable values should be masked'
    )
    assert 'Median telescope pointing offset: -3.0, +2.0 arcsec' in caplog.text, (
        'The measured pointing offset should be reported'
    )


def test_failure_messages(tmp_path, caplog):
    table, rdxdir, truth, coords, flux = write_field(tmp_path, images=STALE_IMAGES, stale=STALE)
    solve_first(rdxdir, truth)
    caplog.set_level(logging.DEBUG, logger='pipnick')
    # In this dense field, the stale images find a few chance matches, which
    # are rejected
    astrometry_all(
        table, rdxdir=rdxdir, solver=FakeSolver(truth, coords),
        catalog=FakeCatalog(coords, flux), wide_search=0., max_submissions=0, qa=False
    )
    assert 'votes' in caplog.text, 'The details of each match should be logged at debug level'
    summary = Table.read(rdxdir / 'astrometry' / 'summary.ecsv')
    unsolved = summary[summary['status'] == 'unsolved']
    assert len(unsolved) == 2 and all('RMS residual' in r for r in unsolved['reason']), (
        'The reason the chance matches were rejected should be recorded'
    )
    assert all(
        'check the pointing' in s and '--max_submissions' in s for s in unsolved['suggestion']
    ), 'The suggestion should cover the rejected matches and why the images were not submitted'

    # With only a few detected sources, too few match
    caplog.clear()
    astrometry_all(
        table, rdxdir=rdxdir, solver=FakeSolver(truth, coords),
        catalog=FakeCatalog(coords, flux), wide_search=0., max_submissions=0, detect_thresh=300.,
        qa=False
    )
    assert 'matched the' in caplog.text and 'minimum is 4' in caplog.text, (
        'Failures to match should give the numbers of sources'
    )
    summary = Table.read(rdxdir / 'astrometry' / 'summary.ecsv')
    unsolved = summary[summary['status'] == 'unsolved']
    assert all('--search' in s for s in unsolved['suggestion']), (
        'Images with too few matches should get a suggestion about the pointing'
    )


def test_qa_and_detection(tmp_path, caplog):
    table, rdxdir, truth, coords, flux = write_field(tmp_path)
    solver = FakeSolver(truth, coords)
    catalog = FakeCatalog(coords, flux)

    # Detection controls and the dry-run QA figures
    astrometry_all(
        table, rdxdir=rdxdir, dry_run=True, qa=True, max_sources=20, detect_thresh=5.
    )
    sources = Table.read(sorted((rdxdir / 'astrometry' / 'dry_run').glob('*.ecsv'))[0])
    assert len(sources) == 20, 'The number of detected sources should be limited'
    assert len(list((rdxdir / 'astrometry' / 'qa').glob('*.png'))) == 4, (
        'A dry run should write a QA figure for every image'
    )
    astrometry_all(table, rdxdir=rdxdir, dry_run=True, detect_thresh=500.)
    assert 'sources were detected in' in caplog.text, 'Images with few sources should be noted'

    # QA figures are redrawn for images solved in the run, and drawn for any
    # image without one
    qa_dir = rdxdir / 'astrometry' / 'qa'
    dry = {f.name: f.stat().st_mtime_ns for f in sorted(qa_dir.glob('*.png'))}
    astrometry_all(table, rdxdir=rdxdir, solver=solver, catalog=catalog)
    solved = {f.name: f.stat().st_mtime_ns for f in sorted(qa_dir.glob('*.png'))}
    assert len(solved) == 4 and all(solved[k] != dry[k] for k in dry), (
        'By default, the QA figure of each image solved in the run should be redrawn'
    )
    assert all(matplotlib.image.imread(f).ndim == 3 for f in sorted(qa_dir.glob('*.png'))), (
        'The QA figures should be readable PNG images'
    )
    missing = sorted(qa_dir.glob('*.png'))[0]
    missing.unlink()
    astrometry_all(table, rdxdir=rdxdir, solver=solver, catalog=catalog)
    rerun = {f.name: f.stat().st_mtime_ns for f in sorted(qa_dir.glob('*.png'))}
    assert missing.name in rerun, 'A missing QA figure should be drawn'
    assert all(rerun[k] == solved[k] for k in rerun if k != missing.name), (
        'The figures of images solved in an earlier run should not be redrawn'
    )
    missing.unlink()
    astrometry_all(table, rdxdir=rdxdir, solver=solver, catalog=catalog, qa=False)
    assert not missing.is_file(), 'No QA figure should be drawn if qa is False'
