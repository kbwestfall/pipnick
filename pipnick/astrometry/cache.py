"""
The cache of the astrometry.net results.

Every submission to astrometry.net is recorded in
``astrometry/astrometry_net.ecsv``, keyed by the hash of what was submitted
(see :func:`anet_input_hash`), and each solution is saved exactly as returned
(``astrometry/anet_<hash>.fits``), so that an image is never submitted twice
and its solution can be restored; see :func:`restore_solution`.
"""

import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

from pipnick import log
from pipnick.astrometry.fit import wcs_from_header
from pipnick.astrometry.header import apply_wcs
from pipnick.astrometry.services import icrs_header
from pipnick.pkg.paths import points_file
from pipnick.pkg.provenance import digest


def anet_input_hash(image_data, sources, hints, upload_image):
    """
    Return a hash of what would be submitted to astrometry.net for an image:
    the source positions (rounded to 0.01 pixel) or, if the image is
    uploaded, the masked image data; and the solver settings.  See
    :func:`~pipnick.pkg.provenance.digest`.
    """
    if upload_image:
        data = np.ascontiguousarray(image_data).tobytes()
    else:
        data = np.round(np.column_stack([sources['x'], sources['y']]), 2).tobytes()
    return digest(data, repr(sorted(hints.items())))



def _anet_index_file(astro_dir):
    """
    Return the index of the astrometry.net submissions; see
    :func:`record_anet`.
    """
    return astro_dir / 'astrometry_net.ecsv'



def _anet_record(astro_dir, input_hash):
    """
    Return the record of the astrometry.net submission of an input, or None;
    see :func:`record_anet`.
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



def record_anet(astro_dir, image, status, submission_id, solution):
    """
    Record the result of an astrometry.net submission.

    The solutions returned by astrometry.net are saved exactly as returned
    (``astrometry/anet_<hash>.fits``, where ``<hash>`` is the hash of the
    input; see :func:`anet_input_hash`), and never modified.  The result of
    each submission is recorded separately, in the index
    ``astrometry/astrometry_net.ecsv``, with one row per input: the input hash,
    the image, the status (``'solved'``, ``'failed'``, ``'rejected'`` by
    :func:`~pipnick.astrometry.solve._check_submission`, or ``'pending'``), the
    submission ID, the file with the solution (``'none'`` if there is none),
    and the date.

    Parameters
    ----------
    astro_dir : :class:`pathlib.Path`
        The directory with the index.
    image : :class:`~pipnick.astrometry.solve.AstrometryImage`
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



def restore_solution(image, images, astro_dir):
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
    image : :class:`~pipnick.astrometry.solve.AstrometryImage`
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

    points = points_file(image.path)
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



def invalidate_chains(images, by_name):
    """
    Invalidate bootstrapped solutions whose reference image is not solved.

    :func:`restore_solution` checks each bootstrapped solution against the
    solution of its reference as it was found on disk.  If the reference is
    itself a bootstrapped solution that was invalidated (and so will be
    solved again), its dependents must be solved again too.  This is applied
    repeatedly, so that the invalidation propagates along chains of
    bootstrapped solutions.

    Parameters
    ----------
    images : list
        The :class:`~pipnick.astrometry.solve.AstrometryImage` objects.  Their states are updated.
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

