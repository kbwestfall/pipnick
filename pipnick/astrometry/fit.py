"""
Fit, construct, and characterize astrometric (WCS) solutions.

All pixel coordinates are 0-based (the ``numpy`` convention, as returned by
``photutils``); the WCS objects follow the FITS (1-based) convention internally,
as usual for :class:`astropy.wcs.WCS`.
"""

import inspect
import warnings

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.wcs import FITSFixedWarning, WCS
from astropy.wcs.utils import fit_wcs_from_points, proj_plane_pixel_scales

from pipnick.pkg.exceptions import PipnickDataModelError


def _fit_wcs_from_points(x, y, coords, sip_degree=None):
    """
    Call :func:`astropy.wcs.utils.fit_wcs_from_points` for 0-based pixel
    coordinates, independent of the ``astropy`` version.

    Versions of ``astropy`` with the ``origin`` argument (added after
    ``astropy`` 8.0) are told that the coordinates are 0-based.  In earlier
    versions, the fit is performed for 0-based coordinates, so the argument is
    not needed.
    """
    if 'origin' in inspect.signature(fit_wcs_from_points).parameters:
        return fit_wcs_from_points((x, y), coords, origin=0, sip_degree=sip_degree)
    return fit_wcs_from_points((x, y), coords, sip_degree=sip_degree)


def wcs_residuals(wcs, x, y, coords):
    """
    Return the angular distance between the positions of a set of sources
    predicted by a WCS and their known coordinates.

    Parameters
    ----------
    wcs : :class:`astropy.wcs.WCS`
        The astrometric solution.
    x, y : :class:`numpy.ndarray`
        The 0-based pixel coordinates of the sources.
    coords : :class:`astropy.coordinates.SkyCoord`
        The celestial coordinates of the sources.

    Returns
    -------
    :class:`numpy.ndarray`
        The residuals in arcsec.
    """
    return wcs.pixel_to_world(x, y).separation(coords).to('arcsec').value


def fit_wcs(x, y, coords, sip_degree=None, sigma=3., min_rejection=0.1, maxiter=5):
    """
    Fit an astrometric solution to sources with known pixel and celestial
    coordinates.

    The solution is a gnomonic (TAN) projection, with optional SIP distortion,
    fit using :func:`astropy.wcs.utils.fit_wcs_from_points`.  Outliers are
    iteratively rejected and the solution refit.

    Parameters
    ----------
    x, y : array-like
        The 0-based pixel coordinates of the sources.
    coords : :class:`astropy.coordinates.SkyCoord`
        The celestial coordinates of the sources.
    sip_degree : int, optional
        The degree of the SIP distortion polynomials.  If None, no distortion
        is fit.
    sigma : float, optional
        Sources with residuals larger than this multiple of the RMS residual
        are rejected.
    min_rejection : float, optional
        Sources with residuals smaller than this (in arcsec) are never
        rejected, so that sources are not rejected from a nearly perfect fit.
    maxiter : int, optional
        Maximum number of rejection iterations.

    Returns
    -------
    wcs : :class:`astropy.wcs.WCS`
        The astrometric solution.
    stats : dict
        The fit statistics: the number of sources used (``npts``), the RMS
        residual of the sources used in arcsec (``rms``), and a boolean array
        selecting the sources used (``used``).

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if there are too few sources to fit the solution: 3 without
        distortion, or :math:`(d+1)(d+2)` for SIP degree :math:`d`.
    """
    _x = np.asarray(x, dtype=float)
    _y = np.asarray(y, dtype=float)
    min_points = 3 if sip_degree is None else (sip_degree + 1) * (sip_degree + 2)
    if len(_x) < min_points:
        raise PipnickDataModelError(
            f'At least {min_points} sources are needed to fit the astrometric solution; '
            f'{len(_x)} provided.'
        )
    used = np.ones(len(_x), dtype=bool)
    for _ in range(maxiter):
        wcs = _fit_wcs_from_points(_x[used], _y[used], coords[used], sip_degree=sip_degree)
        resid = wcs_residuals(wcs, _x, _y, coords)
        rms = float(np.sqrt(np.mean(resid[used]**2)))
        keep = resid <= max(sigma * rms, min_rejection)
        if np.array_equal(keep, used) or np.sum(keep) < min_points:
            break
        used = keep
    return wcs, {'npts': int(np.sum(used)), 'rms': rms, 'used': used}


def scale_and_rotation(wcs):
    """
    Return the pixel scale, rotation, and parity of an astrometric solution.

    Parameters
    ----------
    wcs : :class:`astropy.wcs.WCS`
        The astrometric solution.

    Returns
    -------
    scale : float
        The mean pixel scale in arcsec per pixel.
    rotation : float
        The angle in degrees of the pixel +y axis from north, measured toward
        east for the usual (east-left) parity.
    parity : int
        The sign of the determinant of the pixel scale matrix; -1 for the usual
        orientation of the sky, with east to the left when north is up.
    """
    cd = wcs.celestial.pixel_scale_matrix
    scale = float(np.mean(proj_plane_pixel_scales(wcs.celestial))) * 3600.
    rotation = float(np.degrees(np.arctan2(cd[1,0], cd[1,1])))
    return scale, rotation, int(np.sign(np.linalg.det(cd)))


def predict_wcs(ref_wcs, center, shape):
    """
    Predict the astrometric solution of an image from that of another image
    with the same setup, given the coordinates of its center.

    The prediction has the pixel scale and rotation of the reference
    solution, without distortion, and is centered at the provided
    coordinates.

    Parameters
    ----------
    ref_wcs : :class:`astropy.wcs.WCS`
        The reference solution.
    center : :class:`astropy.coordinates.SkyCoord`
        The coordinates of the center of the image.
    shape : tuple
        The shape of the image (rows, columns).

    Returns
    -------
    :class:`astropy.wcs.WCS`
        The predicted solution.
    """
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.cd = ref_wcs.celestial.pixel_scale_matrix
    wcs.wcs.crval = [center.icrs.ra.deg, center.icrs.dec.deg]
    # FITS pixel coordinates of the image center
    wcs.wcs.crpix = [(shape[1] + 1) / 2, (shape[0] + 1) / 2]
    wcs.pixel_shape = (shape[1], shape[0])
    return wcs


def image_center(wcs, shape):
    """
    Return the celestial coordinates of the center of an image.

    Parameters
    ----------
    wcs : :class:`astropy.wcs.WCS`
        The astrometric solution.
    shape : tuple
        The shape of the image (rows, columns).

    Returns
    -------
    :class:`astropy.coordinates.SkyCoord`
        The coordinates of the image center, always in the ICRS.  Solutions
        found by pipnick are in the ICRS (see
        :func:`~pipnick.astrometry.services.icrs_header`), but a solution
        provided by the user may be in another frame (e.g., FK5), and
        coordinates in different frames cannot be compared directly.
    """
    return to_icrs(wcs.pixel_to_world((shape[1] - 1) / 2, (shape[0] - 1) / 2))


def to_icrs(coords):
    """
    Convert celestial coordinates to the ICRS.

    Unlike :attr:`astropy.coordinates.SkyCoord.icrs`, the returned coordinates
    carry no attributes of the original frame.  Coordinates converted from FK5
    with :attr:`~astropy.coordinates.SkyCoord.icrs` keep the FK5 ``equinox``,
    and :class:`astropy.coordinates.SkyCoord` treats coordinates with
    different attributes as being in different frames (e.g., it refuses to
    compute offsets between them).

    Parameters
    ----------
    coords : :class:`astropy.coordinates.SkyCoord`
        The coordinates.

    Returns
    -------
    :class:`astropy.coordinates.SkyCoord`
        The coordinates in the ICRS.
    """
    icrs = coords.icrs
    return SkyCoord(icrs.ra, icrs.dec, frame='icrs')


def wcs_from_transform(ref_wcs, rotation, shift, shape, ngrid=10):
    r"""
    Construct the astrometric solution of an image from that of a reference
    image, given the rigid transformation between their pixel coordinates.

    The pixel coordinates :math:`{\bf p}` of the image are related to the
    pixel coordinates :math:`{\bf q}` of the reference image by :math:`{\bf
    q} = {\bf R}(\theta) {\bf p} + {\bf t}`, where :math:`{\bf R}(\theta)` is
    the rotation matrix; see :func:`~pipnick.astrometry.match.fit_rigid`.  A
    grid of positions in the image is mapped to the sky through this
    transformation and the reference solution, and a solution is fit to the
    grid; see :func:`fit_wcs`.  Any distortion in the reference solution is
    approximated by a SIP polynomial of the same degree.

    Parameters
    ----------
    ref_wcs : :class:`astropy.wcs.WCS`
        The reference solution.
    rotation : float
        The rotation angle :math:`\theta` in radians.
    shift : array-like
        The shift :math:`{\bf t}` in pixels.
    shape : tuple
        The shape of the image (rows, columns).
    ngrid : int, optional
        The number of grid positions along each axis.

    Returns
    -------
    :class:`astropy.wcs.WCS`
        The solution for the image.
    """
    gx, gy = np.meshgrid(
        np.linspace(0, shape[1] - 1, ngrid), np.linspace(0, shape[0] - 1, ngrid)
    )
    gx = gx.ravel()
    gy = gy.ravel()
    c = np.cos(rotation)
    s = np.sin(rotation)
    qx = c * gx - s * gy + shift[0]
    qy = s * gx + c * gy + shift[1]
    sip_degree = None if ref_wcs.sip is None else int(ref_wcs.sip.a_order)
    wcs, _ = fit_wcs(gx, gy, ref_wcs.pixel_to_world(qx, qy), sip_degree=sip_degree, maxiter=1)
    wcs.pixel_shape = (shape[1], shape[0])
    return wcs


def wcs_from_header(hdr):
    """
    Construct the celestial astrometric solution from a header, or return None
    if it does not have one.

    Parameters
    ----------
    hdr : :class:`astropy.io.fits.Header`
        The header.

    Returns
    -------
    :class:`astropy.wcs.WCS`, None
        The solution, or None if the header does not define a celestial
        solution.
    """
    with warnings.catch_warnings():
        # Headers without data (e.g., cached solutions) give a warning that the
        # WCS has more axes than the image
        warnings.simplefilter('ignore', FITSFixedWarning)
        wcs = WCS(hdr, relax=True)
    if not wcs.has_celestial:
        return None
    return wcs.celestial if wcs.sip is None else wcs

