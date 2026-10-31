"""
Write astrometric solutions to the headers of the reduced images.

A solution replaces any existing WCS in the header (including SIP distortion),
and is recorded with provenance keywords; see :func:`apply_wcs`.
"""

import re
from pathlib import Path

from astropy.io import fits

from pipnick import log
from pipnick.pkg.paths import points_file
from pipnick.pkg.provenance import wcs_hash
from pipnick.pkg.version import version


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
    hdr['PNASTVER'] = (version, 'pipnick version (astrometry)')
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
        hdr['PNASTPOE'] = _pointing_error_card(pointing_error)
    return hdr



def _pointing_error_card(pointing_error):
    """
    Return the value and comment of the ``PNASTPOE`` keyword; see
    :func:`apply_wcs`.
    """
    return round(pointing_error, 3), 'Pointing error in the table (arcmin)'



def set_pointing_error(image, pointing_error):
    """
    Record the pointing error of the astrometric solution of a reduced image.

    Parameters
    ----------
    image : str, :class:`pathlib.Path`
        The reduced image, which is modified in place.
    pointing_error : float, None
        The pointing error in arcmin (``PNASTPOE``; see :func:`apply_wcs`), or
        None to remove it.

    Returns
    -------
    :class:`astropy.io.fits.Header`
        The updated primary header of the image.
    """
    path = Path(image).absolute()
    with fits.open(path, mode='update') as hdu:
        if pointing_error is None:
            hdu[0].header.remove('PNASTPOE', ignore_missing=True)
        else:
            hdu[0].header['PNASTPOE'] = _pointing_error_card(pointing_error)
        hdr = hdu[0].header.copy()
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
        - ``PNASTPOE``: the distance in arcmin between the image center and the
          pointing in the metadata table, if it is larger than the search
          distance of :func:`~pipnick.pipelines.astrometry.astrometry_all`,
          i.e., if the pointing is likely incorrect.

    The solution is written in its own frame (``RADESYS``).  Solutions found by
    :func:`~pipnick.pipelines.astrometry.astrometry_all` are in the ICRS; a
    solution provided by the user should preferably be in the ICRS as well
    (e.g., fit to Gaia positions with :func:`~pipnick.astrometry.fit.fit_wcs`),
    although other frames are handled correctly.

    Solutions provided by the user (``method='points'``) are also saved to
    ``astrometry/<image>_points.fits`` in the directory of the image, so that
    they are restored by :func:`~pipnick.pipelines.astrometry.astrometry_all`
    if the image is rewritten (e.g., by ``pipnick_reduce``).  They are only
    replaced if :func:`~pipnick.pipelines.astrometry.astrometry_all` is run
    with ``resolve=True``.

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
        sidecar = points_file(path)
        sidecar.parent.mkdir(parents=True, exist_ok=True)
        fits.PrimaryHDU(header=solution).writeto(sidecar, overwrite=True)
    log.info(f'Applied the {method} astrometric solution to {path.name}.')
    return hdr

