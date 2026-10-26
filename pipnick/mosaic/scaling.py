r"""
Relative flux scaling of the images of a mosaic.

The images of a mosaic can have different transparency, airmass, and exposure
time, so they are scaled to a common flux scale before they are combined.
The scales are measured from the sources common to pairs of overlapping
images: for each pair :math:`(i, j)` and each common source, the difference of
the instrumental magnitudes is :math:`m_i - m_j = z_i - z_j`, where
:math:`z_i` is the magnitude offset of image :math:`i` with respect to the
reference image (:math:`z_{\rm ref} = 0`).  The offsets are found by a
weighted least-squares fit to all of the differences at once (see
:func:`solve_offsets`), so that images that share no sources with the
reference are still scaled through the images they overlap.  Image :math:`i`
is then multiplied by :math:`10^{0.4 z_i}`.
"""

import numpy as np

from astropy.coordinates import SkyCoord
from scipy import linalg

from pipnick import log


def mutual_matches(ra1, dec1, ra2, dec2, radius):
    """
    Match two lists of positions one-to-one.

    Two positions are matched if each is the nearest to the other, and they
    are separated by less than ``radius``.

    Parameters
    ----------
    ra1, dec1, ra2, dec2 : :class:`numpy.ndarray`
        The ICRS coordinates (deg) of the two lists.
    radius : float
        The match radius (arcsec).

    Returns
    -------
    i1, i2 : :class:`numpy.ndarray`
        The indices of the matched positions in each list.
    """
    if len(ra1) == 0 or len(ra2) == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    c1 = SkyCoord(ra1, dec1, unit='deg')
    c2 = SkyCoord(ra2, dec2, unit='deg')
    near12, sep, _ = c1.match_to_catalog_sky(c2)
    near21, _, _ = c2.match_to_catalog_sky(c1)
    mutual = (near21[near12] == np.arange(len(c1))) & (sep.arcsec < radius)
    return np.where(mutual)[0], near12[mutual]


def _connected(nimages, pairs):
    """
    Return the connected component of each image in the graph whose edges are
    the given pairs of images.
    """
    component = np.arange(nimages)
    for i, j in pairs:
        ci, cj = component[i], component[j]
        if ci != cj:
            component[component == cj] = ci
    return component


def solve_offsets(nimages, measurements, reference=0, floor=0.01, rej=3., maxiter=5, min_common=3):
    r"""
    Solve for the magnitude offsets of a set of images.

    Parameters
    ----------
    nimages : int
        The number of images.
    measurements : list
        For each pair of images, a tuple ``(i, j, dm, err)``, where ``dm`` and
        ``err`` are arrays with the magnitude differences :math:`m_i - m_j` of
        the common sources and their errors.
    reference : int, optional
        The reference image, whose offset is 0.
    floor : float, optional
        A systematic error (mag) added in quadrature to each difference.
    rej : float, optional
        The rejection threshold for the normalized residuals.
    maxiter : int, optional
        The maximum number of rejection iterations.
    min_common : int, optional
        Pairs with fewer common sources are not used.

    Returns
    -------
    offset, error : :class:`numpy.ndarray`
        The offset :math:`z_i` of each image and its error.  Images that are
        not connected to the reference by pairs with at least ``min_common``
        sources have an offset and error of 0.
    scaled : :class:`numpy.ndarray`
        Boolean array flagging the images whose offset was measured (always
        True for the reference).
    nused : :class:`numpy.ndarray`
        The number of differences used for each image, after rejection.
    """
    usable = [m for m in measurements if len(m[2]) >= min_common]
    component = _connected(nimages, [(m[0], m[1]) for m in usable])
    scaled = component == component[reference]
    offset = np.zeros(nimages, dtype=float)
    error = np.zeros(nimages, dtype=float)
    nused = np.zeros(nimages, dtype=int)
    free = [k for k in range(nimages) if scaled[k] and k != reference]
    if len(free) == 0:
        return offset, error, scaled, nused
    column = {k: n for n, k in enumerate(free)}

    # One row per difference
    rows_i = np.concatenate([np.full(len(m[2]), m[0]) for m in usable if scaled[m[0]]])
    rows_j = np.concatenate([np.full(len(m[2]), m[1]) for m in usable if scaled[m[0]]])
    dm = np.concatenate([np.asarray(m[2], dtype=float) for m in usable if scaled[m[0]]])
    err = np.sqrt(
        np.concatenate([np.asarray(m[3], dtype=float) for m in usable if scaled[m[0]]])**2
        + floor**2
    )
    x = np.zeros((len(dm), len(free)), dtype=float)
    for k, n in column.items():
        x[rows_i == k, n] += 1.
        x[rows_j == k, n] -= 1.
    use = np.ones(len(dm), dtype=bool)
    for _ in range(maxiter + 1):
        w = use / err**2
        normal = np.dot(x.T * w, x)
        coeff = linalg.solve(normal, np.dot(x.T * w, dm))
        resid = (dm - np.dot(x, coeff)) / err
        scatter = np.std(resid[use])
        reject = use & (np.abs(resid) > rej * max(scatter, 1.))
        if not np.any(reject):
            break
        use = use & np.logical_not(reject)
    cov = linalg.inv(normal)
    dof = int(np.sum(use)) - len(free)
    if dof > 0:
        rchi2 = np.sum((resid[use])**2) / dof
        if rchi2 > 1.:
            cov *= rchi2
    for k, n in column.items():
        offset[k] = coeff[n]
        error[k] = np.sqrt(cov[n,n])
    for k in range(nimages):
        nused[k] = int(np.sum(use & ((rows_i == k) | (rows_j == k))))
    nrej = int(np.sum(np.logical_not(use)))
    if nrej > 0:
        log.info(f'Rejected {nrej} of {len(dm)} magnitude differences when solving the scales.')
    return offset, error, scaled, nused
