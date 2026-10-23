"""
Astrometric calibration of the reduced science images.

The astrometric solution (WCS) of each image is found with as few submissions
to the astrometry.net web service as possible:

    1. Existing solutions are kept, and solutions are restored from the cache
       of earlier astrometry.net results and from solutions provided by the
       user (see :func:`apply_wcs`), e.g., after the images are rewritten by
       ``pipnick_reduce``.

    2. The images are grouped by their pointing (see
       :func:`group_by_pointing`).  In each group, one image (the reference)
       is solved, and the others are solved by matching their sources to
       those of the reference (*bootstrapping*).

    3. A group without a solved image is solved, if possible, from a solved
       image in another group (if their fields overlap), or by matching its
       sources to the Gaia catalog, using the solution predicted from a solved
       image in the same setup and the telescope pointing offset.

    4. An image that still cannot be solved is matched to the Gaia catalog
       with a *wide search*, which allows for an incorrect telescope pointing
       in the metadata table (e.g., a stale pointing in the raw header), with
       a stricter acceptance threshold.  Images with the same incorrect
       pointing then use the pointing error measured for the first.

    5. Only an image that cannot be solved in any of these ways is submitted
       to astrometry.net.

Solutions are written to the headers of the reduced images, with ``PNAST*``
keywords that record how each was found.

All solutions found by the pipeline are in the ICRS.  astrometry.net labels
its solutions as FK5 (J2000), and they are relabeled as ICRS; see
:func:`~pipnick.astrometry.services.icrs_header` for why this is correct.
"""

import hashlib
import re
import tempfile
from pathlib import Path

import numpy as np

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

from pipnick import log
from pipnick.astrometry.detect import detect_sources
from pipnick.astrometry.fit import (
    fit_wcs, image_center, predict_wcs, scale_and_rotation, to_icrs, wcs_from_header,
    wcs_from_transform, wcs_hash, wcs_residuals
)
from pipnick.astrometry.match import match_frames
from pipnick.astrometry.qa import plot_solution
from pipnick.astrometry.services import (
    AstrometryNetSolver, GaiaCatalog, astrometry_net_hints, icrs_header
)
from pipnick.pipelines.reduction import load_setup, plan_products
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickPathError


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
        :func:`~pipnick.astrometry.detect.detect_sources`.
    hints : dict
        The astrometry.net settings; see
        :func:`~pipnick.astrometry.services.astrometry_net_hints`.
    input_hash : str
        The hash of what would be submitted to astrometry.net, used to find
        cached results; see :func:`_input_hash`.
    comb_id : int, optional
        The combination group of the image.
    """
    def __init__(
        self, path, header, shape, pointing, epoch, scale, fwhm, sources, hints, input_hash,
        comb_id=None
    ):
        self.path = path
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
        # figure; see plot_solution
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
        is applied; see :func:`apply_wcs`.
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
        :func:`_suggestion`).  Failures to solve the image locally are
        ``'few_matches'``, ``'few_matches_wide'``, ``'rejected'``, and
        ``'no_reference'``; failures to submit it are ``'limit'``,
        ``'waiting'``, and ``'submission_failed'``.  The last failure of each
        kind is kept.
    message : str
        The explanation.
    """
    if code in ['limit', 'waiting', 'submission_failed']:
        image.submit_reason = (code, message)
    else:
        image.reason = (code, message)
    log.info(message)


class _Session:
    """
    The settings and the remote services for one run of
    :func:`astrometry_all`, including the count of astrometry.net
    submissions.
    """
    def __init__(
        self, solver, catalog, api_key, astro_dir, upload_image, max_submissions, min_matches,
        max_rms, search, wide_search, wide_min_matches
    ):
        self._solver = solver
        self._api_key = api_key
        self.catalog = catalog
        self.astro_dir = astro_dir
        self.upload_image = upload_image
        self.max_submissions = max_submissions
        self.min_matches = min_matches
        self.max_rms = max_rms
        self.search = search
        self.wide_search = wide_search
        self.wide_min_matches = wide_min_matches
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
        if self.max_submissions is not None and self.nsubmitted >= self.max_submissions:
            _fail(
                image, 'limit',
                f'Not submitting {image.name}: the limit of {self.max_submissions} astrometry.net '
                'submissions has been reached.'
            )
            return False
        return True


def _remove_wcs(hdr):
    """
    Remove the WCS keywords (including any SIP distortion) and the pipnick
    astrometry keywords from a header, in place.
    """
    pattern = re.compile(
        r'^(WCSAXES|CTYPE\d|CUNIT\d|CRVAL\d|CRPIX\d|CDELT\d|CROTA\d|CD\d_\d|PC\d_\d|PV\d_\d+|'
        r'LONPOLE|LATPOLE|RADESYS|RADECSYS|EQUINOX|[AB]P?_ORDER|[AB]P?_\d+_\d+|[AB]_DMAX|'
        r'PNAST\w*)$'
    )
    for key in list(dict.fromkeys(hdr.keys())):
        if pattern.match(key) is not None:
            del hdr[key]


def _solution_header(
    wcs, method, stats=None, reference=None, reference_hash=None, input_hash=None,
    submission_id=None, pointing_error=None
):
    """
    Construct the header cards for an astrometric solution and its provenance;
    see :func:`apply_wcs`.
    """
    hdr = wcs.to_header(relax=True)
    hdr['PNASTMTH'] = (method, 'Astrometric solution method')
    hdr['PNASTWH'] = (wcs_hash(wcs), 'Hash of the astrometric solution')
    if stats is not None:
        hdr['PNASTNPT'] = (stats['npts'], 'Number of sources used in the solution')
        hdr['PNASTRMS'] = (round(stats['rms'], 4), 'RMS residual of the solution (arcsec)')
    if reference is not None:
        hdr['PNASTREF'] = (reference, 'Reference for the solution')
    if reference_hash is not None:
        hdr['PNASTRH'] = (reference_hash, 'Hash of the reference solution')
    if input_hash is not None:
        hdr['PNASTIN'] = (input_hash, 'Hash of the astrometry.net input')
    if submission_id is not None:
        hdr['PNASTSUB'] = (str(submission_id), 'astrometry.net submission ID')
    if pointing_error is not None:
        hdr['PNASTPOE'] = (round(pointing_error, 3), 'Pointing error in the table (arcmin)')
    return hdr


def apply_wcs(
    image, wcs, method='points', stats=None, reference=None, reference_hash=None,
    input_hash=None, submission_id=None, pointing_error=None
):
    """
    Write an astrometric solution to the header of a reduced image.

    Any existing solution is replaced.  The solution is accompanied by
    provenance keywords:

        - ``PNASTMTH``: how the solution was found: ``'astrometry.net'``,
          ``'bootstrap'`` (from another image), ``'gaia'`` (by matching to the
          Gaia catalog), or ``'points'`` (provided by the user, e.g., from
          :func:`~pipnick.astrometry.fit.fit_wcs`);
        - ``PNASTWH``: a hash of the solution;
        - ``PNASTNPT`` and ``PNASTRMS``: the number of sources used and the
          RMS residual in arcsec, if available;
        - ``PNASTREF`` and ``PNASTRH``: the reference image (or catalog) and
          the hash of its solution, for solutions found from another image;
        - ``PNASTIN`` and ``PNASTSUB``: the hash of the astrometry.net input
          and the submission ID, for astrometry.net solutions;
        - ``PNASTPOE``: the distance in arcmin between the image center and
          the pointing in the metadata table, if it is larger than the search
          distance of :func:`astrometry_all`, i.e., if the pointing is likely
          incorrect.

    The solution is written in its own frame (``RADESYS``).  Solutions found
    by :func:`astrometry_all` are in the ICRS; a solution provided by the
    user should preferably be in the ICRS as well (e.g., fit to Gaia
    positions with :func:`~pipnick.astrometry.fit.fit_wcs`), although other
    frames are handled correctly.

    Solutions provided by the user (``method='points'``) are also saved to
    ``astrometry/<image>_points.fits`` in the directory of the image, so that
    they are restored by :func:`astrometry_all` if the image is rewritten
    (e.g., by ``pipnick_reduce``).  They are only replaced if
    :func:`astrometry_all` is run with ``resolve=True``.

    Parameters
    ----------
    image : str, :class:`pathlib.Path`
        The reduced image, which is modified in place.
    wcs : :class:`astropy.wcs.WCS`
        The solution.
    method : str, optional
        How the solution was found; see above.
    stats : dict, optional
        The fit statistics; see :func:`~pipnick.astrometry.fit.fit_wcs`.
    reference : str, optional
        The reference image or catalog.
    reference_hash : str, optional
        The hash of the reference solution.
    input_hash : str, optional
        The hash of the astrometry.net input.
    submission_id : str, optional
        The astrometry.net submission ID.
    pointing_error : float, optional
        The pointing error in arcmin, if the pointing is likely incorrect.

    Returns
    -------
    :class:`astropy.io.fits.Header`
        The updated primary header of the image.
    """
    path = Path(image).absolute()
    solution = _solution_header(
        wcs, method, stats=stats, reference=reference, reference_hash=reference_hash,
        input_hash=input_hash, submission_id=submission_id, pointing_error=pointing_error
    )
    with fits.open(path, mode='update') as hdu:
        _remove_wcs(hdu[0].header)
        hdu[0].header.update(solution)
        hdr = hdu[0].header.copy()
    if method == 'points':
        sidecar = path.parent / 'astrometry' / f'{path.stem}_points.fits'
        sidecar.parent.mkdir(parents=True, exist_ok=True)
        fits.PrimaryHDU(header=solution).writeto(sidecar, overwrite=True)
    log.info(f'Applied the {method} astrometric solution to {path.name}.')
    return hdr


def _input_hash(image_data, sources, hints, upload_image):
    """
    Return a hash of what would be submitted to astrometry.net for an image:
    the source positions (rounded to 0.01 pixel) or, if the image is
    uploaded, the masked image data; and the solver settings.  The hash is the
    first 16 hexadecimal digits of the SHA-1 hash.
    """
    sha = hashlib.sha1()
    if upload_image:
        sha.update(np.ascontiguousarray(image_data).tobytes())
    else:
        sha.update(np.round(np.column_stack([sources['x'], sources['y']]), 2).tobytes())
    sha.update(repr(sorted(hints.items())).encode())
    return sha.hexdigest()[:16]


def _read_image(path):
    """
    Read the data of a reduced image, with masked pixels set to 0, and its mask
    and primary header.
    """
    with fits.open(path, memmap=False) as hdu:
        hdr = hdu[0].header.copy()
        data = np.asarray(hdu[0].data, dtype=float)
        try:
            mask = np.asarray(hdu['MASK'].data, dtype=bool)
        except KeyError:
            mask = np.zeros(data.shape, dtype=bool)
    data[mask] = 0.
    return data, mask, hdr


def _load_images(
    metadata, rdxdir, tag, camera, upload_image, max_sources=100, detect_thresh=5., fwhm=None
):
    """
    Read the reduced science images and detect their sources.

    A warning is issued for images with fewer than 10 detected sources, which
    are likely to be hard to solve.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table.
    rdxdir : :class:`pathlib.Path`
        The reduction directory.
    tag : str
        The setup tag; see :func:`~pipnick.metadata.setup_tag`.
    camera : :class:`~pipnick.cameras.cameras.NickelCamera`
        The camera used to obtain the data.
    upload_image : bool
        Images will be uploaded to astrometry.net, instead of source lists.
    max_sources : int, optional
        The maximum number of (brightest) sources to detect in each image.
    detect_thresh : float, optional
        The detection threshold, as a multiple of the standard deviation of the
        background.
    fwhm : float, optional
        The FWHM of the sources in arcsec.  If None, the typical value for the
        camera is used.

    Returns
    -------
    list
        The :class:`AstrometryImage` objects, ordered by combination group.
        Images that do not exist or have no coordinates are skipped.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if an image is from a different setup.
    """
    images = []
    for product in plan_products(metadata, rdxdir):
        if product['frametype'] != 'science':
            continue
        path = product['file']
        if not path.is_file():
            log.warning(f'{path.name} does not exist; reduce the data first.  Skipping.')
            continue
        data, mask, hdr = _read_image(path)
        if hdr.get('PNSETUP') != tag:
            raise PipnickDataModelError(f'{path.name} is not from setup {tag}.')
        i = product['indx'][0]
        if metadata['ra'][i] is np.ma.masked or metadata['dec'][i] is np.ma.masked:
            log.warning(f'{path.name} has no coordinates.  Skipping.')
            continue
        pointing = SkyCoord(
            metadata['ra'][i], metadata['dec'][i], unit=(units.hourangle, units.deg)
        )
        scale = camera.plate_scale * int(metadata['xbin'][i])
        fwhm_pix = (camera.fwhm_guess if fwhm is None else fwhm) / scale
        sources = detect_sources(
            data, mask=mask, fwhm=fwhm_pix, thresh=detect_thresh, nmax=max_sources
        )
        # Allow for a pointing error of up to the image width
        radius = max(data.shape) * scale / 3600.
        hints = astrometry_net_hints(pointing, radius, scale)
        images.append(
            AstrometryImage(
                path, hdr, data.shape, pointing, Time(metadata['mjd'][i], format='mjd').jyear,
                scale, fwhm_pix, sources, hints, _input_hash(data, sources, hints, upload_image),
                comb_id=product['comb_id']
            )
        )
        log.info(f'{path.name}: detected {len(sources)} sources.')
        if len(sources) < 10:
            log.warning(
                f'Only {len(sources)} sources were detected in {path.name}, which may be hard to '
                'solve; see the troubleshooting advice in the documentation.'
            )
    return images


def _anet_index_file(astro_dir):
    """
    Return the index of the astrometry.net submissions; see
    :func:`_record_anet`.
    """
    return astro_dir / 'astrometry_net.ecsv'


def _anet_record(astro_dir, input_hash):
    """
    Return the record of the astrometry.net submission of an input, or None;
    see :func:`_record_anet`.
    """
    index_file = _anet_index_file(astro_dir)
    if not index_file.is_file():
        return None
    index = Table.read(index_file, format='ascii.ecsv')
    match = index['input_hash'] == input_hash
    if not np.any(match):
        return None
    row = index[np.where(match)[0][0]]
    return {k: str(row[k]) for k in index.colnames}


def _record_anet(astro_dir, image, status, submission_id, solution):
    """
    Record the result of an astrometry.net submission.

    The solutions returned by astrometry.net are saved exactly as returned
    (``astrometry/anet_<hash>.fits``, where ``<hash>`` is the hash of the
    input; see :func:`_input_hash`), and never modified.  The result of each
    submission is recorded separately, in the index
    ``astrometry/astrometry_net.ecsv``, with one row per input: the input
    hash, the image, the status (``'solved'``, ``'failed'``, ``'rejected'``
    by :func:`_check_submission`, or ``'pending'``), the submission ID, the
    file with the solution (``'none'`` if there is none), and the date.

    Parameters
    ----------
    astro_dir : :class:`pathlib.Path`
        The directory with the index.
    image : :class:`AstrometryImage`
        The submitted image.
    status : str
        The status of the submission.
    submission_id : str
        The astrometry.net submission ID.
    solution : str
        The name of the file with the solution, or ``'none'``.
    """
    index_file = _anet_index_file(astro_dir)
    row = {
        'input_hash': image.input_hash, 'image': image.name, 'status': status,
        'submission_id': str(submission_id), 'solution': solution, 'date': Time.now().isot
    }
    if index_file.is_file():
        index = Table.read(index_file, format='ascii.ecsv')
        index.remove_rows(np.where(index['input_hash'] == image.input_hash)[0])
        index.add_row(row)
    else:
        index = Table(rows=[row])
    astro_dir.mkdir(parents=True, exist_ok=True)
    index.write(index_file, format='ascii.ecsv', overwrite=True)


def _restore_solution(image, images, astro_dir):
    """
    Restore the existing solution of an image, if it is still valid, or a
    cached solution.

    An existing solution is kept unless it was bootstrapped from a reference
    image whose solution has since changed, or it is an astrometry.net
    solution not labeled as ICRS (see
    :func:`~pipnick.astrometry.services.icrs_header`).  Otherwise, a solution provided by
    the user (``astrometry/<image>_points.fits``) or a cached astrometry.net
    result (matched by the hash of the input)
    is applied (see :func:`_anet_record`).  A failed or rejected submission
    prevents the image from being submitted again, and a pending submission
    is monitored instead of submitting the image again.

    Parameters
    ----------
    image : :class:`AstrometryImage`
        The image.  Its state is updated.
    images : dict
        All images, keyed by file name.
    astro_dir : :class:`pathlib.Path`
        The directory with the cached solutions.
    """
    method = image.header.get('PNASTMTH')
    if method == 'astrometry.net' and image.header.get('RADESYS') != 'ICRS':
        # Applied before astrometry.net solutions were relabeled as ICRS (see
        # icrs_header); restore it from the cache below, which relabels it
        log.info(f'Relabeling the astrometry.net solution for {image.name} as ICRS.')
        method = None
    if method is not None:
        ref = images.get(image.header.get('PNASTREF'))
        if method == 'bootstrap' and (
            ref is None or ref.solution_hash != image.header.get('PNASTRH')
        ):
            log.info(f'The reference solution for {image.name} has changed.')
            return
        image.set_solution(image.header)
        log.info(f'{image.name} already has a {method} solution.')
        return

    points = astro_dir / f'{image.path.stem}_points.fits'
    if points.is_file():
        image.set_solution(apply_wcs(image.path, wcs_from_header(fits.getheader(points))))
        return
    record = _anet_record(astro_dir, image.input_hash)
    if record is None:
        return
    if record['status'] == 'pending':
        log.info(f'The astrometry.net submission for {image.name} is pending.')
        image.status = 'pending'
        image.submission_id = record['submission_id']
    elif record['status'] == 'solved':
        log.info(f'Using the cached astrometry.net solution for {image.name}.')
        # The solution as returned by astrometry.net, relabeled as ICRS
        solution = fits.Header.fromfile(astro_dir / record['solution'])
        image.set_solution(
            apply_wcs(
                image.path, wcs_from_header(icrs_header(solution)), method='astrometry.net',
                input_hash=image.input_hash, submission_id=record['submission_id']
            )
        )
    else:
        log.info(
            f'An earlier astrometry.net submission of {image.name} failed; it will not be '
            'submitted again unless resolve is set.'
        )
        image.submit_failed = True


def _invalidate_chains(images, by_name):
    """
    Invalidate bootstrapped solutions whose reference image is not solved.

    :func:`_restore_solution` checks each bootstrapped solution against the
    solution of its reference as it was found on disk.  If the reference is
    itself a bootstrapped solution that was invalidated (and so will be
    solved again), its dependents must be solved again too.  This is applied
    repeatedly, so that the invalidation propagates along chains of
    bootstrapped solutions.

    Parameters
    ----------
    images : list
        The :class:`AstrometryImage` objects.  Their states are updated.
    by_name : dict
        The same images, keyed by file name.
    """
    changed = True
    while changed:
        changed = False
        for image in images:
            if image.status != 'solved' or image.method != 'bootstrap':
                continue
            ref = by_name.get(image.header.get('PNASTREF'))
            if ref is None or ref.status != 'solved':
                log.info(f'The reference image for {image.name} will be solved again.')
                image.status = None
                image.wcs = None
                image.method = None
                changed = True


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


def _pointing_offset(images):
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
    the pointing offset of the solved images; see :func:`_pointing_offset`.

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
    dlon, dlat = _pointing_offset(images)
    return image.pointing.spherical_offsets_by(dlon, dlat)


def _pointing_error(image, wcs, session):
    """
    Check the pointing in the metadata table against a solution.

    If the center of the solved image is farther from its pointing in the
    metadata table than ``search``, the pointing is likely incorrect (e.g., a
    stale pointing in the raw header).  This is reported in a warning, and the
    image is excluded from the pointing offset used to predict the positions
    of other images (see :func:`_pointing_offset`).

    Returns
    -------
    float, None
        The distance in arcmin between the image center and the pointing, if
        it is larger than ``search``, or None otherwise.
    """
    center = image_center(wcs, image.shape)
    offset = center.separation(image.pointing).arcsec
    if offset <= session.search:
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
    _min_matches = session.min_matches if min_matches is None else min_matches
    scale, rotation, parity = scale_and_rotation(wcs)
    ref_scale, ref_rotation, ref_parity = scale_and_rotation(ref_wcs)
    drot = (rotation - ref_rotation + 180.) % 360. - 180.
    reasons = []
    if stats['npts'] < _min_matches:
        reasons.append(f'{stats["npts"]} matched sources (minimum is {_min_matches})')
    if stats['rms'] > session.max_rms:
        reasons.append(f'RMS residual is {stats["rms"]:.2f} arcsec (maximum is {session.max_rms})')
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
    :func:`~pipnick.astrometry.qa.plot_solution`.
    """
    xref, yref = wcs.world_to_pixel(coords)
    return {
        'x': np.asarray(image.sources['x'])[match['indx']],
        'y': np.asarray(image.sources['y'])[match['indx']], 'xref': xref, 'yref': yref,
        'label': label
    }


def _log_pointing_offset(images):
    """
    Log the median telescope pointing offset (see :func:`_pointing_offset`),
    which sets the predicted positions of the images.
    """
    n = len([i for i in images if i.status == 'solved' and i.pointing_error is None])
    if n == 0:
        return
    dlon, dlat = _pointing_offset(images)
    log.info(
        f'Median telescope pointing offset: {dlon.value:+.1f}, {dlat.value:+.1f} arcsec (east, '
        f'north) from {n} solved images.'
    )


def _suggestion(image, images):
    """
    Suggest what the user can do to solve an image that could not be solved,
    based on why it could not be solved locally and why it could not be
    submitted; see :func:`_fail`.
    """
    suggestions = {
        'few_matches': 'check the pointing in the metadata table, or increase --search',
        'few_matches_wide': (
            'check the pointing in the metadata table, or increase --search or --wide_search'
        ),
        'rejected': (
            'the matches may be by chance if the pointing is wrong, so check the pointing in '
            'the metadata table; otherwise, consider --max_rms'
        ),
        'no_reference': (
            'solve another image in its pointing group or setup first'
            if any(i.status == 'solved' for i in images)
            else 'no image in the setup is solved, so a submission is needed'
        ),
        'limit': 'allow a submission (--max_submissions)',
        'waiting': 'rerun to check the pending submission',
        'submission_failed': (
            'try --upload_image or provide a solution; use --resolve and --images to resubmit'
        )
    }
    found = []
    if len(image.sources) < 10:
        found.append(
            'few detected sources: try --detect_thresh, --upload_image, or provide a solution '
            '(fit_wcs and apply_wcs)'
        )
    for reason in [image.reason, image.submit_reason]:
        if reason is not None:
            found.append(suggestions[reason[0]])
    return 'see the log' if len(found) == 0 else '; '.join(found)


def _summarize(images, astro_dir):
    """
    Log a summary of the solution of each image, and write it to
    ``astrometry/summary.ecsv``.

    The summary gives, for each image, its status (``solved``, ``unsolved``,
    ``pending``, or ``skipped`` if not selected), the method, the reference,
    the number of matched sources, the RMS residual (arcsec), and the
    pointing error (arcmin; see :func:`_pointing_error`) of its solution, and,
    for images that are not solved, the reason and a suggestion.
    """
    rows = []
    for image in images:
        status = (
            image.status if image.status in ['solved', 'pending', 'skipped'] else 'unsolved'
        )
        hdr = image.header if image.status == 'solved' else {}
        reason = '' if status == 'solved' else ' '.join(
            r[1] for r in [image.reason, image.submit_reason] if r is not None
        )
        suggestion = _suggestion(image, images) if status == 'unsolved' else ''
        rows.append(
            {
                'image': image.name, 'status': status, 'method': hdr.get('PNASTMTH', '-'),
                'reference': hdr.get('PNASTREF', '-'), 'nmatch': hdr.get('PNASTNPT', -1),
                'rms': hdr.get('PNASTRMS', -1.), 'pointing_error': hdr.get('PNASTPOE', -1.),
                'reason': reason if len(reason) > 0 else '-',
                'suggestion': suggestion if len(suggestion) > 0 else '-'
            }
        )
    tbl = Table(rows=rows)
    # Missing values are masked, not set to an invalid number
    for name in ['nmatch', 'rms', 'pointing_error']:
        tbl[name] = np.ma.MaskedArray(tbl[name], mask=np.asarray(tbl[name]) < 0)
    astro_dir.mkdir(parents=True, exist_ok=True)
    tbl.write(astro_dir / 'summary.ecsv', format='ascii.ecsv', overwrite=True)
    log.info(
        'Astrometry summary (also in astrometry/summary.ecsv):\n'
        + '\n'.join(tbl['image', 'status', 'method', 'reference', 'nmatch', 'rms',
                         'pointing_error'].pformat(max_lines=-1, max_width=-1))
    )
    for row in tbl:
        if row['status'] == 'unsolved':
            log.warning(f'{row["image"]} is not solved: {row["reason"]}  Suggestion: '
                        f'{row["suggestion"]}.')


def _write_qa(images, astro_dir, previous=None):
    """
    Write a QA figure for images to ``astrometry/qa/<image>.png``; see
    :func:`~pipnick.astrometry.qa.plot_solution`.

    Parameters
    ----------
    images : list
        The :class:`AstrometryImage` objects.
    astro_dir : :class:`pathlib.Path`
        The ``astrometry`` directory.
    previous : set, optional
        The names of the images that were solved before this run, which are
        skipped, along with images that are not solved.  If None, a figure is
        written for every image (e.g., for a dry run).
    """
    qa_dir = astro_dir / 'qa'
    qa_dir.mkdir(parents=True, exist_ok=True)
    for image in images:
        if previous is not None and (image.status != 'solved' or image.name in previous):
            continue
        data, mask, _ = _read_image(image.path)
        title = f'{image.name}: {len(image.sources)} detected sources'
        if image.status == 'solved':
            title += f'; {image.method} solution'
        plot_solution(
            data, mask, image.sources, qa_dir / f'{image.path.stem}.png', matched=image.matched,
            title=title
        )
    log.info(f'QA figures written to {qa_dir}.')


def _select(images, selection):
    """
    Select the images to solve.

    Parameters
    ----------
    images : list
        The :class:`AstrometryImage` objects.  Their ``selected`` attribute is
        set.
    selection : list
        The images to select, by file name, file name without the extension,
        or combination group (``comb_id``).

    Raises
    ------
    ValueError
        Raised if an entry does not match any image.
    """
    unknown = []
    for image in images:
        image.selected = False
    for item in selection:
        found = False
        for image in images:
            if (
                str(item) in [image.name, image.path.stem]
                or (str(item).isdigit() and image.comb_id == int(item))
            ):
                image.selected = True
                found = True
        if not found:
            unknown.append(str(item))
    if len(unknown) > 0:
        raise ValueError(f'No reduced science image matches: {", ".join(unknown)}.')


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
        search=session.search / image.scale, tol=max(2., image.fwhm)
    )
    _log_match(image, ref.name, match)
    nmatch = 0 if match is None else len(match['indx'])
    if nmatch < session.min_matches:
        _fail(
            image, 'few_matches',
            f'{nmatch} of the {len(image.sources)} sources in {image.name} matched the '
            f'{len(ref.sources)} sources in {ref.name} (minimum is {session.min_matches}).'
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
    search = (session.wide_search * 60. if wide else session.search) / image.scale
    min_matches = (
        max(session.min_matches, session.wide_min_matches) if wide else session.min_matches
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


def _handle_result(image, status, header, submission_id, session):
    """
    Handle the result of an astrometry.net submission: save and apply a
    solution, or record a failure or a pending submission; see
    :func:`_record_anet`.

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
        _record_anet(session.astro_dir, image, 'pending', submission_id, 'none')
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
        _record_anet(session.astro_dir, image, result, submission_id, solution)
        image.submit_failed = True
        return
    _record_anet(session.astro_dir, image, 'solved', submission_id, solution)
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
    if session.upload_image:
        log.info(f'Uploading {image.name} to astrometry.net.')
        data, _, hdr = _read_image(image.path)
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
    _handle_result(image, status, header, submission_id, session)


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
    if session.wide_search > 0:
        return _bootstrap_gaia(image, ref, images, session, wide=True)
    return False


def _group_reference(group):
    """
    Return the solved image to use as the reference for a pointing group,
    preferring solutions that were not bootstrapped, or None.
    """
    solved = [i for i in group if i.status == 'solved']
    if len(solved) == 0:
        return None
    direct = [i for i in solved if i.method != 'bootstrap']
    return direct[0] if len(direct) > 0 else solved[0]


def _solve_group(group, images, session):
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
    ref = _group_reference(group)
    if ref is None:
        if any(i.status == 'pending' for i in group):
            # Wait for the pending submission
            return
        candidates = [i for i in group if i.status is None]
        if len(candidates) == 0:
            return
        candidates.sort(key=_nsources, reverse=True)
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


def _nsources(image):
    """
    Return the number of detected sources in an image.
    """
    return len(image.sources)


def _report_dry_run(images, groups, astro_dir):
    """
    Report what :func:`astrometry_all` would do, and write the source lists
    that would be submitted to astrometry.net.
    """
    outdir = astro_dir / 'dry_run'
    outdir.mkdir(parents=True, exist_ok=True)
    lines = ['Dry run; nothing is submitted.  Pointing groups:']
    nunsolved = 0
    for i, group in enumerate(groups):
        lines.append(f'    Group {i+1}:')
        for image in group:
            state = (
                f'solved ({image.method})' if image.status == 'solved'
                else 'not selected' if image.status == 'skipped'
                else 'not submitted (earlier submission failed)' if image.submit_failed
                else 'unsolved'
            )
            lines.append(f'        {image.name}: {len(image.sources)} sources; {state}')
            tbl = Table(
                {'x': image.sources['x'], 'y': image.sources['y'], 'flux': image.sources['flux']},
                meta={'hints': image.hints, 'input_hash': image.input_hash}
            )
            tbl.write(outdir / f'{image.path.stem}_sources.ecsv', overwrite=True)
        if _group_reference(group) is None:
            nunsolved += 1
            candidates = sorted(group, key=_nsources, reverse=True)
            lines.append(f'        Reference if a submission is needed: {candidates[0].name}')
    order = [
        max(g, key=_nsources).name for g in sorted(groups, key=_group_order)
        if _group_reference(g) is None
    ]
    if len(order) > 0:
        lines.append(f'    Order of any submissions: {", ".join(order)}')
    nsolved = len([i for i in images if i.status == 'solved'])
    minimum = 1 if nsolved == 0 and nunsolved > 0 else 0
    lines.append(
        f'{nunsolved} groups without a solved image.  At least {minimum} and, if all '
        f'bootstrapping fails, at most {len([i for i in images if i.status is None])} '
        'astrometry.net submissions would be made.'
    )
    lines.append(f'Source lists written to {outdir}.')
    log.info('\n'.join(lines))


def astrometry_all(
    table, rdxdir=None, solver=None, catalog=None, api_key=None, resolve=False,
    upload_image=False, dry_run=False, max_submissions=None, group_radius=1., min_matches=4,
    max_rms=0.3, search=60., wide_search=10., wide_min_matches=10, images=None, qa=False,
    max_sources=100, detect_thresh=5., fwhm=None
):
    """
    Find the astrometric solutions of the reduced science images.

    See the module docstring for the strategy, which minimizes the number of
    submissions to astrometry.net.  The solutions are written to the image
    headers; see :func:`apply_wcs`.  The cached astrometry.net results, Gaia
    queries, and user-provided solutions are kept in the ``astrometry``
    subdirectory of the reduction directory.

    Parameters
    ----------
    table : str, :class:`pathlib.Path`
        The per-setup metadata file used to reduce the images; see
        :func:`~pipnick.pipelines.reduction.reduce_all`.
    rdxdir : str, :class:`pathlib.Path`, optional
        The reduction directory.  If None, the directory saved in the metadata
        file is used.
    solver : object, optional
        The astrometry.net solver.  If None,
        :class:`~pipnick.astrometry.services.AstrometryNetSolver` is
        constructed when it is first needed.
    catalog : object, optional
        The source catalog.  If None,
        :class:`~pipnick.astrometry.services.GaiaCatalog` is used, with its
        queries cached in the ``astrometry`` directory.
    api_key : str, optional
        The astrometry.net API key, if ``solver`` is None; see
        :class:`~pipnick.astrometry.services.AstrometryNetSolver`.
    resolve : bool, optional
        Ignore the existing and cached solutions, and solve the images again
        (only the selected images, if ``images`` is provided).  This can
        submit jobs to astrometry.net for images that were solved before.
    upload_image : bool, optional
        Upload images to astrometry.net, instead of the lists of their
        detected sources.
    dry_run : bool, optional
        Detect the sources, group the images, and report what would be
        submitted, without solving or submitting anything.
    max_submissions : int, optional
        The maximum number of submissions to astrometry.net.  If None, there
        is no limit.
    group_radius : float, optional
        The maximum distance in arcmin between the pointings of images in a
        group; see :func:`group_by_pointing`.
    min_matches : int, optional
        The minimum number of matched sources for a bootstrapped solution.
    max_rms : float, optional
        The maximum RMS residual in arcsec for a bootstrapped solution.
    search : float, optional
        The maximum error in arcsec in the expected offset between matched
        images (or between an image and the catalog), e.g., from the error in
        the telescope pointing.
    wide_search : float, optional
        For images that cannot be solved otherwise, the maximum error in
        arcmin in the expected position of the image when matching to the Gaia
        catalog, to allow for an incorrect pointing in the metadata table.  If
        0, this wide search is not done.
    wide_min_matches : int, optional
        The minimum number of matched sources for a solution found by the wide
        search.
    images : list, optional
        Only solve these images, given by file name, file name without the
        extension, or combination group (``comb_id``).  The other images are
        not changed, but their solutions are used as references.  If None, all
        images are solved.
    qa : bool, optional
        Write a QA figure for each image solved in this run (for every image,
        for a dry run) to ``astrometry/qa``; see
        :func:`~pipnick.astrometry.qa.plot_solution`.
    max_sources : int, optional
        The maximum number of (brightest) sources detected in each image, and
        submitted to astrometry.net.
    detect_thresh : float, optional
        The source detection threshold, as a multiple of the standard
        deviation of the background.
    fwhm : float, optional
        The FWHM of the sources in arcsec, for their detection.  If None, the
        typical value for the camera is used.

    Changing ``max_sources``, ``detect_thresh``, or ``fwhm`` changes the
    sources that would be submitted, so cached astrometry.net results for the
    previous sources are not reused.

    A summary of the solution of each image is logged and written to
    ``astrometry/summary.ecsv``; see :func:`_summarize`.

    Returns
    -------
    list
        The paths to the solved images.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickPathError`
        Raised if the reduction directory does not exist.
    """
    metafile, metadata, _rdxdir, tag, camera = load_setup(table, rdxdir=rdxdir)
    if not _rdxdir.is_dir():
        raise PipnickPathError(f'{_rdxdir} does not exist; reduce the data first.')
    astro_dir = _rdxdir / 'astrometry'
    log.info(f'Astrometric calibration of the images in {_rdxdir}')

    _images = _load_images(
        metadata, _rdxdir, tag, camera, upload_image, max_sources=max_sources,
        detect_thresh=detect_thresh, fwhm=fwhm
    )
    if len(_images) == 0:
        log.warning('No reduced science images to calibrate.')
        return []
    if images is not None:
        _select(_images, images)
    by_name = {image.name: image for image in _images}
    for image in _images:
        if not resolve or not image.selected:
            _restore_solution(image, by_name, astro_dir)
    _invalidate_chains(_images, by_name)
    for image in _images:
        if not image.selected and image.status != 'solved':
            image.status = 'skipped'
    _log_pointing_offset(_images)
    groups = group_by_pointing(_images, group_radius)

    if dry_run:
        _report_dry_run(_images, groups, astro_dir)
        if qa:
            _write_qa(_images, astro_dir)
        return [image.path for image in _images if image.status == 'solved']

    session = _Session(
        solver, GaiaCatalog(astro_dir) if catalog is None else catalog, api_key, astro_dir,
        upload_image, max_submissions, min_matches, max_rms, search, wide_search, wide_min_matches
    )
    session.pending = set(image.name for image in _images if image.status == 'pending')
    previous = set(image.name for image in _images if image.status == 'solved')
    # Solve the groups with a solved image first, then the richest groups
    groups.sort(key=_group_order)
    for group in groups:
        _solve_group(group, _images, session)

    # Check the pending submissions, and then solve the remaining groups
    for image in _images:
        if image.status != 'pending':
            continue
        status, header, submission_id = session.solver().monitor(image.submission_id)
        _handle_result(image, status, header, submission_id, session)
        if image.status == 'pending':
            log.warning(
                f'The astrometry.net submission {image.submission_id} for {image.name} has not '
                'finished; rerun to check it again.'
            )
    groups.sort(key=_group_order)
    for group in groups:
        _solve_group(group, _images, session)

    if qa:
        _write_qa(_images, astro_dir, previous=previous)
    _log_pointing_offset(_images)
    _summarize(_images, astro_dir)
    solved = [image for image in _images if image.status == 'solved']
    methods = {}
    for image in solved:
        methods[image.method] = methods.get(image.method, 0) + 1
    log.info(
        f'Solved {len(solved)} of {len(_images)} images ('
        + ', '.join(f'{n} {m}' for m, n in sorted(methods.items()))
        + f'), with {session.nsubmitted} astrometry.net submissions.'
    )
    return [image.path for image in solved]


def _group_order(group):
    """
    Return the sort key used to order the pointing groups: groups with a
    solved image first, and then by decreasing number of sources in their
    richest image, so that any submission is of the image most likely to be
    solved.
    """
    return (0 if _group_reference(group) is not None else 1, -max(_nsources(i) for i in group))
