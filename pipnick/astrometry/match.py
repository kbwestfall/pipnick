"""
Match sources between two sets of pixel coordinates, related by a rigid
transformation (a shift and a small rotation, with no change of scale).

All pixel coordinates are 0-based.
"""

import numpy as np
from scipy.spatial import cKDTree


def apply_rigid(x, y, rotation, shift):
    r"""
    Apply a rigid transformation to a set of pixel coordinates.

    The transformation is :math:`{\bf q} = {\bf R}(\theta) {\bf p} + {\bf t}`,
    where :math:`{\bf R}(\theta)` is the rotation matrix.

    Parameters
    ----------
    x, y : :class:`numpy.ndarray`
        The pixel coordinates, :math:`{\bf p}`.
    rotation : float
        The rotation angle, :math:`\theta`, in radians.
    shift : array-like
        The shift, :math:`{\bf t}`, in pixels.

    Returns
    -------
    tuple
        The transformed coordinates, :math:`{\bf q}`.
    """
    c = np.cos(rotation)
    s = np.sin(rotation)
    return c * x - s * y + shift[0], s * x + c * y + shift[1]


def fit_rigid(x, y, xref, yref):
    r"""
    Fit the rigid transformation between two sets of matched pixel
    coordinates.

    Finds the rotation, :math:`\theta`, and shift, :math:`{\bf t}`, that
    minimize the squared distance between :math:`{\bf R}(\theta) {\bf p} +
    {\bf t}` and :math:`{\bf q}`; see :func:`apply_rigid`.

    Parameters
    ----------
    x, y : :class:`numpy.ndarray`
        The pixel coordinates to transform, :math:`{\bf p}`.
    xref, yref : :class:`numpy.ndarray`
        The matched reference coordinates, :math:`{\bf q}`.

    Returns
    -------
    rotation : float
        The rotation angle in radians.
    shift : :class:`numpy.ndarray`
        The shift in pixels.
    """
    px = x - np.mean(x)
    py = y - np.mean(y)
    qx = xref - np.mean(xref)
    qy = yref - np.mean(yref)
    rotation = float(np.arctan2(np.sum(px * qy - py * qx), np.sum(px * qx + py * qy)))
    mx, my = apply_rigid(np.mean(x), np.mean(y), rotation, [0., 0.])
    return rotation, np.array([np.mean(xref) - mx, np.mean(yref) - my])


def vote_shift(x, y, xref, yref, expected=(0., 0.), search=50., binsize=2.):
    """
    Find the shift between two sets of pixel coordinates by voting.

    The shift between every pair of positions, one from each set, is
    computed, and the shift shared by the most pairs is found from the peak of
    their 2D histogram.  This is robust to sources in one set that are not in
    the other, but assumes any rotation between the sets is small.

    Parameters
    ----------
    x, y : :class:`numpy.ndarray`
        The pixel coordinates to shift.
    xref, yref : :class:`numpy.ndarray`
        The reference pixel coordinates.
    expected : array-like, optional
        The expected shift, from the reference to the shifted coordinates.
    search : float, optional
        Only shifts within this distance (in pixels, along each axis) of the
        expected shift are considered.
    binsize : float, optional
        The histogram bin size in pixels, which should be comparable to the
        uncertainty in the positions.

    Returns
    -------
    shift : :class:`numpy.ndarray`, None
        The shift, such that ``xref ~ x + shift[0]``, or None if no pairs are
        within the search region.
    nvotes : int
        The number of pairs consistent with the shift.
    """
    dx = (xref[None,:] - x[:,None]).ravel() - expected[0]
    dy = (yref[None,:] - y[:,None]).ravel() - expected[1]
    indx = (np.absolute(dx) <= search) & (np.absolute(dy) <= search)
    if not np.any(indx):
        return None, 0
    nbin = max(int(np.ceil(2 * search / binsize)), 1)
    hist, xedges, yedges = np.histogram2d(
        dx[indx], dy[indx], bins=nbin, range=[[-search, search], [-search, search]]
    )
    i, j = np.unravel_index(np.argmax(hist), hist.shape)
    cx = (xedges[i] + xedges[i+1]) / 2
    cy = (yedges[j] + yedges[j+1]) / 2
    near = indx & (np.absolute(dx - cx) <= binsize) & (np.absolute(dy - cy) <= binsize)
    shift = np.array([np.mean(dx[near]) + expected[0], np.mean(dy[near]) + expected[1]])
    return shift, int(np.sum(near))


def match_nearest(x, y, xref, yref, tol):
    """
    Match each position to the nearest reference position, one-to-one.

    Parameters
    ----------
    x, y : :class:`numpy.ndarray`
        The pixel coordinates to match, already transformed to the reference
        frame.
    xref, yref : :class:`numpy.ndarray`
        The reference pixel coordinates.
    tol : float
        The maximum distance in pixels for a match.

    Returns
    -------
    indx : :class:`numpy.ndarray`
        The indices of the matched positions.
    indx_ref : :class:`numpy.ndarray`
        The indices of the matched reference positions.  If more than one
        position is closest to the same reference position, only the closest
        is kept.
    """
    if len(x) == 0 or len(xref) == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    dist, near = cKDTree(np.column_stack([xref, yref])).query(np.column_stack([x, y]))
    indx = np.where(dist <= tol)[0]
    # Keep only the closest match to each reference position
    srt = indx[np.argsort(dist[indx])]
    _, first = np.unique(near[srt], return_index=True)
    indx = np.sort(srt[first])
    return indx, near[indx]


def match_frames(
    x, y, xref, yref, expected=(0., 0.), search=50., tol=3., sigma=3., maxiter=5
):
    """
    Match two sets of sources related by a rigid transformation, and fit the
    transformation.

    The shift is first found by voting (:func:`vote_shift`), and the sources
    are then iteratively matched (:func:`match_nearest`) and the rigid
    transformation fit (:func:`fit_rigid`), rejecting outliers.

    Parameters
    ----------
    x, y : :class:`numpy.ndarray`
        The pixel coordinates of the sources to match.
    xref, yref : :class:`numpy.ndarray`
        The pixel coordinates of the reference sources.
    expected : array-like, optional
        The expected shift; see :func:`vote_shift`.
    search : float, optional
        The search distance for the shift; see :func:`vote_shift`.
    tol : float, optional
        The maximum distance in pixels for a match, after applying the
        transformation.
    sigma : float, optional
        Matches with residuals larger than this multiple of the RMS residual
        (and larger than half of ``tol``) are rejected.
    maxiter : int, optional
        Maximum number of iterations.

    Returns
    -------
    dict, None
        The ``rotation`` (radians) and ``shift`` (pixels) of the
        transformation (see :func:`apply_rigid`), the indices of the matched
        sources (``indx`` and ``indx_ref``), the RMS residual of the matches in
        pixels (``rms``), and the number of pairs that voted for the initial
        shift (``nvotes``; see :func:`vote_shift`).  None if fewer than 2
        sources are matched.
    """
    shift, nvotes = vote_shift(x, y, xref, yref, expected=expected, search=search, binsize=tol)
    if shift is None:
        return None
    rotation = 0.
    for _ in range(maxiter):
        tx, ty = apply_rigid(x, y, rotation, shift)
        indx, indx_ref = match_nearest(tx, ty, xref, yref, tol)
        if len(indx) < 2:
            return None
        rotation, shift = fit_rigid(x[indx], y[indx], xref[indx_ref], yref[indx_ref])
        tx, ty = apply_rigid(x[indx], y[indx], rotation, shift)
        resid = np.hypot(tx - xref[indx_ref], ty - yref[indx_ref])
        rms = float(np.sqrt(np.mean(resid**2)))
        keep = resid <= max(sigma * rms, tol / 2)
        if np.all(keep):
            break
        if np.sum(keep) < 2:
            return None
        rotation, shift = fit_rigid(
            x[indx][keep], y[indx][keep], xref[indx_ref][keep], yref[indx_ref][keep]
        )
    return {
        'rotation': rotation, 'shift': shift, 'indx': indx, 'indx_ref': indx_ref, 'rms': rms,
        'nvotes': nvotes
    }
