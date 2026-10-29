"""
The strategy for finding the astrometric solutions of a set of images.

The images are grouped by pointing (see :func:`group_by_pointing`), and the
groups are solved in order (see :func:`solve_group`): from an image already
solved, by matching their sources to those of a solved image or to the Gaia
catalog, and only otherwise by submitting an image to astrometry.net (see
:func:`handle_result`); see :mod:`pipnick.pipelines.astrometry`.
"""

import tempfile
from pathlib import Path

import numpy as np

from astropy import units
from astropy.io import fits

from pipnick import log
from pipnick.astrometry.cache import record_anet
from pipnick.astrometry.fit import (
    fit_wcs, image_center, predict_wcs, scale_and_rotation, to_icrs, wcs_from_header,
    wcs_from_transform, wcs_residuals
)
from pipnick.astrometry.header import apply_wcs
from pipnick.astrometry.match import match_frames
from pipnick.astrometry.services import AstrometryNetSolver, icrs_header
from pipnick.image.io import read_reduced


def read_image(path):
    """
    Read the data of a reduced image, with masked pixels set to 0, and its mask
    and primary header; see :func:`~pipnick.image.io.read_reduced`.

    The masked pixels are set to 0 because the image uploaded to
    astrometry.net (with ``upload_image``), and the sources detected in it,
    are part of the key of the cached astrometry.net results (see
    :func:`~pipnick.astrometry.cache.anet_input_hash`), so they must not
    change.
    """
    data, mask, _, hdr = read_reduced(path)
    data[mask] = 0.
    return data, mask, hdr



class AstrometryImage:
    """
    The state of one reduced science image during the astrometric
    calibration.

    Parameters
    ----------
    path : :class:`pathlib.Path`
        The reduced image.
    header : :class:`astropy.io.fits.Header`
        Its primary header.
    shape : tuple
        The shape of the image (rows, columns).
    pointing : :class:`astropy.coordinates.SkyCoord`
        The telescope pointing, from the metadata table.
    epoch : float
        The epoch of the observation (Julian year).
    scale : float
        The nominal pixel scale in arcsec per pixel.
    fwhm : float
        The nominal FWHM of the point sources in pixels.
    sources : :class:`astropy.table.Table`
        The detected sources; see
        :func:`~pipnick.image.detect.find_sources`.
    hints : dict
        The astrometry.net settings; see
        :func:`~pipnick.astrometry.services.astrometry_net_hints`.
    input_hash : str
        The hash of what would be submitted to astrometry.net, used to find
        cached results; see :func:`~pipnick.astrometry.cache.anet_input_hash`.
    comb_id : int, optional
        The combination group of the image.
    image : :class:`~pipnick.datamodel.image.ReducedImage`, optional
        The reduced image, with its product and data set.
    """
    def __init__(
        self, path, header, shape, pointing, epoch, scale, fwhm, sources, hints, input_hash,
        comb_id=None, image=None
    ):
        self.path = path
        self.image = image
        self.comb_id = comb_id
        self.header = header
        self.shape = shape
        self.pointing = pointing
        self.epoch = epoch
        self.scale = scale
        self.fwhm = fwhm
        self.sources = sources
        self.hints = hints
        self.input_hash = input_hash
        # The solution state
        self.status = None
        self.wcs = None
        self.method = None
        self.submission_id = None
        self.submit_failed = False
        # The distance (arcmin) between the image center and the pointing in
        # the metadata table, if larger than the search distance (i.e., the
        # pointing is likely incorrect); see _pointing_error
        self.pointing_error = None
        # Whether the image is selected for solving
        self.selected = True
        # The last reasons the image could not be solved locally and could not
        # be submitted, each as a (code, message) tuple; see _fail
        self.reason = None
        self.submit_reason = None
        # The matched sources of a solution found in this run, for the QA
        # figure; see plot_sources
        self.matched = None

    @property
    def name(self):
        """The file name of the image."""
        return self.path.name

    @property
    def solution_hash(self):
        """The hash of the applied solution, or None."""
        return self.header.get('PNASTWH')

    def set_solution(self, header):
        """
        Set the solution state from the header of the image after a solution
        is applied; see :func:`~pipnick.astrometry.header.apply_wcs`.
        """
        self.header = header
        self.wcs = wcs_from_header(header)
        self.method = header['PNASTMTH']
        self.pointing_error = header.get('PNASTPOE')
        self.status = 'solved'
        self.reason = None
        self.submit_reason = None



def _fail(image, code, message):
    """
    Record and log why an image could not be solved.

    Parameters
    ----------
    image : :class:`AstrometryImage`
        The image.
    code : str
        The kind of failure, used to suggest what the user can do (see
        :func:`~pipnick.pipelines.astrometry._suggestion`).  Failures to solve
        the image locally are ``'few_matches'``, ``'few_matches_wide'``,
        ``'rejected'``, and ``'no_reference'``; failures to submit it are
        ``'limit'``, ``'waiting'``, and ``'submission_failed'``.  The last
        failure of each kind is kept.
    message : str
        The explanation.
    """
    if code in ['limit', 'waiting', 'submission_failed']:
        image.submit_reason = (code, message)
    else:
        image.reason = (code, message)
    log.info(message)



class Session:
    """
    The settings and the remote services for one run of
    :func:`~pipnick.pipelines.astrometry.astrometry_all`, including the count of
    astrometry.net submissions.

    Parameters
    ----------
    solver : object
        The astrometry.net solver, or None to construct the default solver
        when it is first needed; see :meth:`solver`.
    catalog : object
        The source catalog, e.g., a
        :class:`~pipnick.astrometry.services.GaiaCatalog`.
    api_key : str
        The astrometry.net API key, for the default solver.
    astro_dir : :class:`pathlib.Path`
        The ``astrometry`` directory.
    settings : :class:`~pipnick.pkg.settings.AstrometrySettings`
        The settings.
    """
    def __init__(self, solver, catalog, api_key, astro_dir, settings):
        self._solver = solver
        self._api_key = api_key
        self.catalog = catalog
        self.astro_dir = astro_dir
        self.settings = settings
        self.nsubmitted = 0
        # The images with pending submissions
        self.pending = set()

    def solver(self):
        """
        Return the astrometry.net solver, constructing the default solver
        only when it is first needed.
        """
        if self._solver is None:
            self._solver = AstrometryNetSolver(api_key=self._api_key)
        return self._solver

    def can_submit(self, image):
        """
        Check if an image can be submitted to astrometry.net.
        """
        if image.status == 'pending':
            return False
        if image.submit_failed:
            _fail(
                image, 'submission_failed',
                f'Not submitting {image.name}: an earlier astrometry.net submission failed.'
            )
            return False
        if len(self.pending) > 0:
            # A pending solution may allow this image to be solved without a
            # submission
            _fail(
                image, 'waiting',
                f'Not yet submitting {image.name}: waiting for the pending astrometry.net '
                'submissions.'
            )
            return False
        limit = self.settings.max_submissions
        if limit is not None and self.nsubmitted >= limit:
            _fail(
                image, 'limit',
                f'Not submitting {image.name}: the limit of {limit} astrometry.net submissions '
                'has been reached.'
            )
            return False
        return True



def group_by_pointing(images, radius):
    """
    Group images by their telescope pointing.

    Each image joins the first group whose first image is within ``radius``
    of its pointing, or starts a new group.

    Parameters
    ----------
    images : list
        The :class:`AstrometryImage` objects.
    radius : float
        The maximum distance in arcmin between pointings in a group.

    Returns
    -------
    list
        The groups, each a list of images.
    """
    groups = []
    for image in images:
        for group in groups:
            if group[0].pointing.separation(image.pointing).arcmin <= radius:
                group.append(image)
                break
        else:
            groups.append([image])
    return groups



def pointing_offset(images):
    """
    Return the median offset between the solved and the telescope pointing of
    the solved images, as longitude and latitude offsets in arcsec; see
    :func:`astropy.coordinates.SkyCoord.spherical_offsets_to`.  Images with
    likely incorrect pointings (see :func:`_pointing_error`) are excluded.
    The offset is 0 if no image is solved.
    """
    dlon = []
    dlat = []
    for image in images:
        if image.status != 'solved' or image.pointing_error is not None:
            continue
        lon, lat = image.pointing.spherical_offsets_to(image_center(image.wcs, image.shape))
        dlon.append(lon.to('arcsec').value)
        dlat.append(lat.to('arcsec').value)
    if len(dlon) == 0:
        return 0. * units.arcsec, 0. * units.arcsec
    return np.median(dlon) * units.arcsec, np.median(dlat) * units.arcsec



def _predicted_center(image, images):
    """
    Predict the coordinates of the center of an image from its pointing and
    the pointing offset of the solved images; see :func:`pointing_offset`.

    If a solved image with a likely incorrect pointing in the metadata table
    (see :func:`_pointing_error`) has the same pointing (to within 1 arcsec),
    the pointing is likely stale for both, and the offset measured for that
    image is used instead.
    """
    for other in images:
        if (
            other.status == 'solved' and other.pointing_error is not None
            and other.pointing.separation(image.pointing).arcsec < 1.
        ):
            dlon, dlat = other.pointing.spherical_offsets_to(image_center(other.wcs, other.shape))
            return image.pointing.spherical_offsets_by(dlon, dlat)
    dlon, dlat = pointing_offset(images)
    return image.pointing.spherical_offsets_by(dlon, dlat)



def _pointing_error(image, wcs, session):
    """
    Check the pointing in the metadata table against a solution.

    If the center of the solved image is farther from its pointing in the
    metadata table than ``search``, the pointing is likely incorrect (e.g., a
    stale pointing in the raw header).  This is reported in a warning, and the
    image is excluded from the pointing offset used to predict the positions
    of other images (see :func:`pointing_offset`).

    Returns
    -------
    float, None
        The distance in arcmin between the image center and the pointing, if
        it is larger than ``search``, or None otherwise.
    """
    center = image_center(wcs, image.shape)
    offset = center.separation(image.pointing).arcsec
    if offset <= session.settings.search:
        return None
    log.warning(
        f'{image.name} is centered at {center.to_string("hmsdms", sep=":", precision=1)}, '
        f'{offset/60:.1f} arcmin from its pointing in the metadata table.  The pointing in the '
        'raw header may be incorrect; consider correcting the RA and Dec in the metadata table '
        'and rerunning pipnick_reduce.'
    )
    return offset / 60



def _accept(image, wcs, stats, ref_wcs, session, min_matches=None):
    """
    Check a solution found from a reference solution (another image or the
    prediction for the Gaia match).

    The solution must use at least ``min_matches`` sources (the session value
    if None), have an RMS residual of at most ``max_rms``, and agree with the
    reference in pixel scale (within 2%), rotation (within 1 degree), and
    parity.
    """
    _min_matches = session.settings.min_matches if min_matches is None else min_matches
    scale, rotation, parity = scale_and_rotation(wcs)
    ref_scale, ref_rotation, ref_parity = scale_and_rotation(ref_wcs)
    drot = (rotation - ref_rotation + 180.) % 360. - 180.
    reasons = []
    if stats['npts'] < _min_matches:
        reasons.append(f'{stats["npts"]} matched sources (minimum is {_min_matches})')
    if stats['rms'] > session.settings.max_rms:
        reasons.append(
            f'RMS residual is {stats["rms"]:.2f} arcsec (maximum is {session.settings.max_rms})'
        )
    if abs(scale / ref_scale - 1) > 0.02 or abs(drot) > 1. or parity != ref_parity:
        reasons.append('the pixel scale or rotation differs from the reference')
    if len(reasons) > 0:
        _fail(image, 'rejected', f'Rejected the solution for {image.name}: {"; ".join(reasons)}.')
        return False
    return True



def _log_match(image, reference, match):
    """
    Log the details of an attempt to match the sources in an image to a
    reference, at the debug level.
    """
    if match is None:
        log.debug(f'{image.name} vs {reference}: no shift found in the search region.')
        return
    log.debug(
        f'{image.name} vs {reference}: shift ({match["shift"][0]:.1f}, {match["shift"][1]:.1f}) '
        f'pixels, rotation {np.degrees(match["rotation"]):.3f} deg, {match["nvotes"]} votes, '
        f'{len(match["indx"])} matches, RMS {match["rms"]:.2f} pixels.'
    )



def _matched_positions(image, match, wcs, coords, label):
    """
    Return the matched sources of a solution for the QA figure: their pixel
    coordinates, the pixel coordinates of their matches predicted by the
    solution, and a label for the reference; see
    :func:`~pipnick.qa.sources.plot_sources`.
    """
    xref, yref = wcs.world_to_pixel(coords)
    return {
        'x': np.asarray(image.sources['x'])[match['indx']],
        'y': np.asarray(image.sources['y'])[match['indx']], 'xref': xref, 'yref': yref,
        'label': label
    }



def log_pointing_offset(images):
    """
    Log the median telescope pointing offset (see :func:`pointing_offset`),
    which sets the predicted positions of the images.
    """
    n = len([i for i in images if i.status == 'solved' and i.pointing_error is None])
    if n == 0:
        return
    dlon, dlat = pointing_offset(images)
    log.info(
        f'Median telescope pointing offset: {dlon.value:+.1f}, {dlat.value:+.1f} arcsec (east, '
        f'north) from {n} solved images.'
    )



def _fit_matched(image, match, ref_wcs, coords):
    """
    Fit the solution of an image from its sources matched to reference
    sources with known coordinates.

    With at least 6 matches, the solution is fit directly to the matched
    sources (see :func:`~pipnick.astrometry.fit.fit_wcs`); otherwise, the
    solution is the reference solution transformed by the fitted shift and
    rotation (see :func:`~pipnick.astrometry.fit.wcs_from_transform`), which
    keeps the reference pixel scale.

    Returns
    -------
    wcs : :class:`astropy.wcs.WCS`
        The solution.
    stats : dict
        The number of sources used and their RMS residual in arcsec.
    """
    x = np.asarray(image.sources['x'])[match['indx']]
    y = np.asarray(image.sources['y'])[match['indx']]
    if len(x) >= 6:
        wcs, stats = fit_wcs(x, y, coords)
        return wcs, {'npts': stats['npts'], 'rms': stats['rms']}
    wcs = wcs_from_transform(ref_wcs, match['rotation'], match['shift'], image.shape)
    resid = wcs_residuals(wcs, x, y, coords)
    return wcs, {'npts': len(x), 'rms': float(np.sqrt(np.mean(resid**2)))}



def _bootstrap_frame(image, ref, images, session):
    """
    Solve an image by matching its sources to those of a solved reference
    image.

    The expected shift between the images is set by their pointings.  The
    solution is accepted subject to :func:`_accept`.

    Returns
    -------
    bool
        True if the image was solved.
    """
    center = _predicted_center(image, images)
    px, py = ref.wcs.world_to_pixel(center)
    expected = (float(px) - (image.shape[1] - 1) / 2, float(py) - (image.shape[0] - 1) / 2)
    match = match_frames(
        np.asarray(image.sources['x']), np.asarray(image.sources['y']),
        np.asarray(ref.sources['x']), np.asarray(ref.sources['y']), expected=expected,
        search=session.settings.search / image.scale, tol=max(2., image.fwhm)
    )
    _log_match(image, ref.name, match)
    nmatch = 0 if match is None else len(match['indx'])
    if nmatch < session.settings.min_matches:
        _fail(
            image, 'few_matches',
            f'{nmatch} of the {len(image.sources)} sources in {image.name} matched the '
            f'{len(ref.sources)} sources in {ref.name} (minimum is '
            f'{session.settings.min_matches}).'
        )
        return False
    # In the ICRS, even if the reference solution is in another frame
    coords = to_icrs(
        ref.wcs.pixel_to_world(
            np.asarray(ref.sources['x'])[match['indx_ref']],
            np.asarray(ref.sources['y'])[match['indx_ref']]
        )
    )
    wcs, stats = _fit_matched(image, match, ref.wcs, coords)
    if not _accept(image, wcs, stats, ref.wcs, session):
        return False
    image.set_solution(
        apply_wcs(
            image.path, wcs, method='bootstrap', stats=stats, reference=ref.name,
            reference_hash=ref.solution_hash, pointing_error=_pointing_error(image, wcs, session)
        )
    )
    image.matched = _matched_positions(image, match, wcs, coords, ref.name)
    log.info(
        f'Solved {image.name} from {ref.name} with {stats["npts"]} sources (RMS '
        f'{stats["rms"]:.2f} arcsec).'
    )
    return True



def _bootstrap_gaia(image, ref, images, session, wide=False):
    """
    Solve an image by matching its sources to the Gaia catalog, starting from
    the solution predicted from a solved image in the same setup (see
    :func:`~pipnick.astrometry.fit.predict_wcs`).

    The solution is accepted subject to :func:`_accept`.

    If ``wide`` is True, the match allows for an error of up to
    ``wide_search`` arcmin in the predicted position of the image, e.g.,
    because the pointing in the metadata table is incorrect, and at least
    ``wide_min_matches`` matched sources are required.

    For any solution, a pointing error larger than ``search`` is recorded in
    the header (``PNASTPOE``) and reported; see :func:`_pointing_error`.

    Returns
    -------
    bool
        True if the image was solved.
    """
    if session.catalog is None:
        return False
    settings = session.settings
    search = (settings.wide_search * 60. if wide else settings.search) / image.scale
    min_matches = (
        max(settings.min_matches, settings.wide_min_matches) if wide else settings.min_matches
    )
    predicted = predict_wcs(ref.wcs, _predicted_center(image, images), image.shape)
    radius = 0.5 * np.hypot(*image.shape) * image.scale / 3600. + search * image.scale / 3600.
    coords, _ = session.catalog.sources(
        image_center(predicted, image.shape), radius, image.epoch
    )
    gx, gy = predicted.world_to_pixel(coords)
    # Catalog sources that could be in the image.  The wide search must allow
    # for the large error in the predicted position; otherwise, sources outside
    # the predicted image are excluded, because in sparse fields they only add
    # false matches.
    margin = search if wide else 0.
    inside = (
        (gx >= -margin) & (gx <= image.shape[1] - 1 + margin) & (gy >= -margin)
        & (gy <= image.shape[0] - 1 + margin)
    )
    match = match_frames(
        np.asarray(image.sources['x']), np.asarray(image.sources['y']), gx[inside], gy[inside],
        search=search, tol=max(2., image.fwhm)
    )
    _log_match(image, 'Gaia DR3' + (' (wide search)' if wide else ''), match)
    nmatch = 0 if match is None else len(match['indx'])
    if nmatch < min_matches:
        _fail(
            image, 'few_matches_wide' if wide else 'few_matches',
            f'{nmatch} of the {len(image.sources)} sources in {image.name} matched the '
            f'{np.sum(inside)} Gaia sources expected in or near the field'
            + (' in the wide search' if wide else '') + f' (minimum is {min_matches}).'
        )
        return False
    wcs, stats = _fit_matched(image, match, predicted, coords[inside][match['indx_ref']])
    if not _accept(image, wcs, stats, predicted, session, min_matches=min_matches):
        return False
    image.set_solution(
        apply_wcs(
            image.path, wcs, method='gaia', stats=stats, reference='Gaia DR3',
            pointing_error=_pointing_error(image, wcs, session)
        )
    )
    image.matched = _matched_positions(
        image, match, wcs, coords[inside][match['indx_ref']], 'Gaia DR3'
    )
    log.info(
        f'Solved {image.name} using Gaia with {stats["npts"]} sources (RMS '
        f'{stats["rms"]:.2f} arcsec).'
    )
    return True



def _check_submission(image, wcs):
    """
    Check that an astrometry.net solution is consistent with the submitted
    settings: its pixel scale must be within the submitted limits, and its
    center within the submitted radius of the telescope pointing.
    """
    scale, _, _ = scale_and_rotation(wcs)
    offset = image_center(wcs, image.shape).separation(image.pointing).deg
    if scale < image.hints['scale_lower'] or scale > image.hints['scale_upper']:
        log.warning(f'The astrometry.net solution for {image.name} has the wrong pixel scale.')
        return False
    if offset > image.hints['radius']:
        log.warning(f'The astrometry.net solution for {image.name} is too far from the pointing.')
        return False
    return True



def handle_result(image, status, header, submission_id, session):
    """
    Handle the result of an astrometry.net submission: save and apply a
    solution, or record a failure or a pending submission; see
    :func:`~pipnick.astrometry.cache.record_anet`.

    A returned solution is saved exactly as returned.  The solution applied to
    the image is a copy relabeled as ICRS (see
    :func:`~pipnick.astrometry.services.icrs_header`), after checking that it
    is consistent with the submitted settings (see :func:`_check_submission`).
    """
    if status == 'pending':
        log.info(
            f'The astrometry.net solution for {image.name} is not yet finished; it will be '
            'checked again later.'
        )
        image.status = 'pending'
        image.submission_id = submission_id
        session.pending.add(image.name)
        record_anet(session.astro_dir, image, 'pending', submission_id, 'none')
        return
    image.status = None
    session.pending.discard(image.name)
    wcs = None
    solution = 'none'
    if status == 'solved':
        solution = f'anet_{image.input_hash}.fits'
        session.astro_dir.mkdir(parents=True, exist_ok=True)
        # The solution exactly as returned by astrometry.net
        header.tofile(session.astro_dir / solution, overwrite=True)
        # astrometry.net labels its solutions as FK5, but they are ICRS; see icrs_header
        wcs = wcs_from_header(icrs_header(header))
        if not _check_submission(image, wcs):
            wcs = None
    if wcs is None:
        log.warning(f'astrometry.net did not solve {image.name}.')
        image.submit_reason = ('submission_failed', f'astrometry.net did not solve {image.name}.')
        result = 'rejected' if status == 'solved' else 'failed'
        record_anet(session.astro_dir, image, result, submission_id, solution)
        image.submit_failed = True
        return
    record_anet(session.astro_dir, image, 'solved', submission_id, solution)
    image.set_solution(
        apply_wcs(
            image.path, wcs, method='astrometry.net', input_hash=image.input_hash,
            submission_id=submission_id, pointing_error=_pointing_error(image, wcs, session)
        )
    )



def _submit(image, session):
    """
    Submit an image to astrometry.net, either as a list of its sources or, if
    ``upload_image`` is set, as an image with its masked pixels set to 0.
    """
    solver = session.solver()
    session.nsubmitted += 1
    if session.settings.upload_image:
        log.info(f'Uploading {image.name} to astrometry.net.')
        data, _, hdr = read_image(image.path)
        with tempfile.TemporaryDirectory() as tmp:
            upload = Path(tmp) / image.name
            fits.PrimaryHDU(data=data, header=hdr).writeto(upload)
            status, header, submission_id = solver.solve_image(upload, image.hints)
    else:
        log.info(f'Submitting the {len(image.sources)} sources in {image.name} to astrometry.net.')
        status, header, submission_id = solver.solve_sources(
            np.asarray(image.sources['x']), np.asarray(image.sources['y']), image.shape[1],
            image.shape[0], image.hints
        )
    handle_result(image, status, header, submission_id, session)



def _nearest_solved(image, images):
    """
    Return the best solved image to use as the reference for another image,
    or None.

    Among the solved images that overlap substantially with the image (with
    centers within half the image width of its predicted center; see
    :func:`_predicted_center`), the nearest one whose solution was not itself
    bootstrapped is preferred, to avoid accumulating errors along chains of
    bootstrapped solutions.  Otherwise, the solved image nearest the predicted
    center is returned.
    """
    solved = [i for i in images if i.status == 'solved' and i is not image]
    if len(solved) == 0:
        return None
    center = _predicted_center(image, images)
    sep = np.array([image_center(i.wcs, i.shape).separation(center).deg for i in solved])
    width = max(image.shape) * image.scale / 3600.
    direct = np.array([i.method != 'bootstrap' for i in solved]) & (sep < width / 2)
    if np.any(direct):
        indx = np.where(direct)[0]
        return solved[int(indx[np.argmin(sep[indx])])]
    return solved[int(np.argmin(sep))]



def _solve_without_reference(image, images, session):
    """
    Solve an image without a solved image in its pointing group: from the
    nearest solved image if their fields overlap, by matching to the Gaia
    catalog, or, failing that, by matching to the Gaia catalog with the wide
    search (see :func:`_bootstrap_gaia`).

    Returns
    -------
    bool
        True if the image was solved.
    """
    ref = _nearest_solved(image, images)
    if ref is None:
        _fail(image, 'no_reference', f'No image in the setup is solved to match {image.name} to.')
        return False
    width = max(image.shape) * image.scale / 3600.
    center = _predicted_center(image, images)
    if image_center(ref.wcs, ref.shape).separation(center).deg < width / 2:
        if _bootstrap_frame(image, ref, images, session):
            return True
    if _bootstrap_gaia(image, ref, images, session):
        return True
    if session.settings.wide_search > 0:
        return _bootstrap_gaia(image, ref, images, session, wide=True)
    return False



def group_reference(group):
    """
    Return the solved image to use as the reference for a pointing group,
    preferring solutions that were not bootstrapped, or None.
    """
    solved = [i for i in group if i.status == 'solved']
    if len(solved) == 0:
        return None
    direct = [i for i in solved if i.method != 'bootstrap']
    return direct[0] if len(direct) > 0 else solved[0]



def solve_group(group, images, session):
    """
    Solve the images in a pointing group.

    If no image in the group is solved, the image with the most detected
    sources is solved from other groups (see
    :func:`_solve_without_reference`) or, failing that, submitted to
    astrometry.net.  The other images are then bootstrapped (see
    :func:`_bootstrap_frame`) from the best overlapping solved image (see
    :func:`_nearest_solved`) or from the reference, falling back to the same
    sequence.
    """
    ref = group_reference(group)
    if ref is None:
        if any(i.status == 'pending' for i in group):
            # Wait for the pending submission
            return
        candidates = [i for i in group if i.status is None]
        if len(candidates) == 0:
            return
        candidates.sort(key=nsources, reverse=True)
        best = candidates[0]
        if not _solve_without_reference(best, images, session) and session.can_submit(best):
            _submit(best, session)
        if best.status != 'solved':
            for image in candidates[1:]:
                image.reason = (
                    'no_reference', f'No image in its pointing group was solved ({best.name} '
                    'was tried first).'
                )
            return
        ref = best
    for image in group:
        if image.status is not None:
            continue
        # Bootstrap from the best overlapping solved image (see
        # _nearest_solved), which need not be the group reference (e.g., if the
        # pointing in the metadata table is incorrect)
        nearest = _nearest_solved(image, images)
        if nearest is not None and nearest is not ref:
            if _bootstrap_frame(image, nearest, images, session):
                continue
        if _bootstrap_frame(image, ref, images, session):
            continue
        if _solve_without_reference(image, images, session):
            continue
        if session.can_submit(image):
            _submit(image, session)



def nsources(image):
    """
    Return the number of detected sources in an image.
    """
    return len(image.sources)



def group_order(group):
    """
    Return the sort key used to order the pointing groups: groups with a
    solved image first, and then by decreasing number of sources in their
    richest image, so that any submission is of the image most likely to be
    solved.
    """
    return (0 if group_reference(group) is not None else 1, -max(nsources(i) for i in group))

