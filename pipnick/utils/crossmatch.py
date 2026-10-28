"""
Cross-match sources by their celestial coordinates.
"""

import numpy as np

from astropy.coordinates import SkyCoord

from pipnick.pkg.exceptions import PipnickParameterError


def match_sky(coords1, coords2, radius, unique='mutual'):
    """
    Match the positions in one list to those in another.

    Each position in ``coords1`` is matched to its nearest neighbor in
    ``coords2``, if they are separated by less than ``radius``.  The
    ``unique`` option sets how positions in ``coords2`` that are the nearest
    neighbor of more than one position in ``coords1`` are treated.

    Parameters
    ----------
    coords1, coords2 : :class:`astropy.coordinates.SkyCoord`
        The positions to match.
    radius : float
        The match radius (arcsec).
    unique : str, optional
        - ``'mutual'``: a match is kept only if each position is the nearest
          neighbor of the other (one-to-one);
        - ``'second'``: a position in ``coords2`` is kept only in its nearest
          match (one-to-one);
        - ``'none'``: all matches are kept, so several positions in
          ``coords1`` can be matched to the same position in ``coords2``.

    Returns
    -------
    i1, i2 : :class:`numpy.ndarray`
        The indices of the matched positions in each list, in the order of
        ``i1``.
    sep : :class:`numpy.ndarray`
        The separations of the matches (arcsec).

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if ``unique`` is not one of the above.
    """
    if unique not in ('mutual', 'second', 'none'):
        raise PipnickParameterError(f'Unknown matching option: {unique}.')
    if len(coords1) == 0 or len(coords2) == 0:
        return np.array([], dtype=int), np.array([], dtype=int), np.array([])
    near, sep, _ = coords1.match_to_catalog_sky(coords2)
    sep = sep.arcsec
    if unique == 'mutual':
        back, _, _ = coords2.match_to_catalog_sky(coords1)
        i1 = np.where((back[near] == np.arange(len(coords1))) & (sep < radius))[0]
        return i1, near[i1], sep[i1]
    i1 = np.where(sep < radius)[0]
    i2 = near[i1]
    sep = sep[i1]
    if unique == 'second':
        # Keep only the nearest match of each position in coords2
        order = np.argsort(sep)
        _, first = np.unique(i2[order], return_index=True)
        keep = np.sort(order[first])
        return i1[keep], i2[keep], sep[keep]
    return i1, i2, sep


def link_sources(ra, dec, radius):
    """
    Identify the same sources in several lists of positions.

    The lists are processed in order.  The positions in each list are matched
    to the sources identified in the previous lists, within ``radius``, and
    take the identifier of the nearest one; positions that do not match are
    added as new sources.  A position is matched to the first position of
    each source, and several positions in one list can be matched to the
    same source.

    Parameters
    ----------
    ra, dec : list
        The right ascension and declination (deg) of the positions in each
        list, as :class:`numpy.ma.MaskedArray` objects.  Masked positions are
        neither matched nor added as new sources; e.g., a list without
        celestial coordinates is fully masked.
    radius : float
        The match radius (arcsec).

    Returns
    -------
    list
        The source identifier of each position in each list, as
        :class:`numpy.ndarray` objects.  Unmatched (and masked) positions get
        new, consecutive identifiers.
    """
    ref_id = []
    ref_ra = []
    ref_dec = []
    ids = []
    next_id = 0
    for _ra, _dec in zip(ra, dec):
        _ra = np.ma.asarray(_ra, dtype=float)
        _dec = np.ma.asarray(_dec, dtype=float)
        sid = np.arange(next_id, next_id + len(_ra))
        next_id += len(_ra)
        ids.append(sid)
        good = np.where(np.logical_not(np.ma.getmaskarray(_ra) | np.ma.getmaskarray(_dec)))[0]
        new = good
        if len(ref_id) > 0 and len(good) > 0:
            ref = SkyCoord(ref_ra, ref_dec, unit='deg')
            coo = SkyCoord(_ra[good].filled(0.), _dec[good].filled(0.), unit='deg')
            i1, i2, _ = match_sky(coo, ref, radius, unique='none')
            sid[good[i1]] = np.asarray(ref_id)[i2]
            new = np.delete(good, i1)
        ref_id += [int(s) for s in sid[new]]
        ref_ra += [float(r) for r in _ra[new]]
        ref_dec += [float(d) for d in _dec[new]]
    return ids
