"""
The Landolt UBVRI photometric standard stars.

The catalog distributed with pipnick (``pipnick/data/catalogs/landolt.ecsv``)
merges two catalogs of Johnson-Kron-Cousins UBVRI standard stars:

    - Landolt (2009, AJ, 137, 4186), VizieR catalog ``J/AJ/137/4186``: 595
      stars around the celestial equator, re-observing the stars of Landolt
      (1992, AJ, 104, 340), with UCAC2 coordinates and proper motions;
    - Landolt (2013, AJ, 146, 131), VizieR catalog ``J/AJ/146/131``: 349 stars
      at northern declinations and the celestial equator, with UCAC4
      coordinates and proper motions.

The catalog is constructed by :func:`build_catalog`, which retrieves both
catalogs from VizieR.  It is distributed with pipnick, so that using it never
requires network access.

This research has made use of the VizieR catalogue access tool, CDS,
Strasbourg, France (Ochsenbein et al. 2000, A&AS, 143, 23).
"""

import importlib.resources
from pathlib import Path

import numpy as np

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.table import MaskedColumn, Table, vstack
from astroquery.vizier import Vizier

from pipnick import log
from pipnick.astrometry.services import propagate_proper_motion
from pipnick.pkg.exceptions import PipnickParameterError
from pipnick.utils.crossmatch import match_sky


def default_colors():
    """
    Return the default color used to calibrate each band: U-B for U, B-V for
    B and V, V-R for R, and V-I for I.
    """
    return {'U': 'U-B', 'B': 'B-V', 'V': 'B-V', 'R': 'V-R', 'I': 'V-I'}


def landolt_band(filtername):
    """
    Return the Landolt band of a filter.

    Parameters
    ----------
    filtername : str
        The filter name.

    Returns
    -------
    str
        The band, or None if the filter is not one of UBVRI.
    """
    return filtername if filtername in ['U', 'B', 'V', 'R', 'I'] else None


def catalog_path():
    """
    Return the path to the Landolt catalog distributed with pipnick.
    """
    return Path(str(importlib.resources.files('pipnick') / 'data' / 'catalogs' / 'landolt.ecsv'))


def _ingest(tbl, names, source):
    """
    Construct the standardized table of one Landolt catalog.

    Parameters
    ----------
    tbl : :class:`astropy.table.Table`
        The catalog as returned by VizieR.
    names : dict
        The VizieR column names of the star name (``name``), the V magnitude
        (``V``) and the colors (``B-V``, ``U-B``, ``V-R``, ``R-I``, ``V-I``),
        each with its error (prefixed by ``e_``), and the numbers of
        observations (``nobs``) and nights (``nnight``).
    source : str
        The catalog reference, recorded for each star.

    Returns
    -------
    :class:`astropy.table.Table`
        The standardized catalog; see :func:`build_catalog`.
    """
    coo = SkyCoord(tbl['RAJ2000'], tbl['DEJ2000'], unit=(units.hourangle, units.deg))
    out = Table()
    out['name'] = [str(n).strip() for n in tbl[names['name']]]
    out['ra'] = MaskedColumn(coo.ra.deg, unit='deg')
    out['dec'] = MaskedColumn(coo.dec.deg, unit='deg')
    # NOTE: The UCAC proper motions in right ascension include the cos(dec)
    # factor, as expected by propagate_proper_motion
    out['pmra'] = MaskedColumn(np.ma.asarray(tbl['pmRA'], dtype=float), unit='mas/yr')
    out['pmdec'] = MaskedColumn(np.ma.asarray(tbl['pmDE'], dtype=float), unit='mas/yr')
    out['ref_epoch'] = np.full(len(tbl), 2000.)
    for key in ['V', 'B-V', 'U-B', 'V-R', 'R-I', 'V-I']:
        out[key] = MaskedColumn(np.ma.asarray(tbl[names[key]], dtype=float), unit='mag')
        out[f'e_{key}'] = MaskedColumn(
            np.ma.asarray(tbl[names[f'e_{key}']], dtype=float), unit='mag'
        )
    out['nobs'] = np.asarray(tbl[names['nobs']], dtype=int)
    out['nnight'] = np.asarray(tbl[names['nnight']], dtype=int)
    out['source'] = [source] * len(tbl)
    return out


def derive_magnitudes(tbl):
    r"""
    Add the UBVRI magnitudes derived from the V magnitude and the colors.

    The magnitudes are :math:`B = V + (B-V)`, :math:`U = B + (U-B)`,
    :math:`R = V - (V-R)`, and :math:`I = V - (V-I)`.  Their errors are
    propagated naively, i.e., ignoring the correlations between the colors.
    A magnitude is masked if any of the values it is derived from are masked.

    Parameters
    ----------
    tbl : :class:`astropy.table.Table`
        The catalog, with the columns ``V``, ``B-V``, ``U-B``, ``V-R``,
        ``V-I`` and their errors (``e_V``, etc.).  The columns ``U``, ``B``,
        ``R``, ``I`` and their errors are added in place, and ``V`` is
        unchanged.
    """
    col = {
        name: np.ma.asarray(tbl[name], dtype=float)
        for name in ['V', 'B-V', 'U-B', 'V-R', 'V-I', 'e_V', 'e_B-V', 'e_U-B', 'e_V-R', 'e_V-I']
    }
    derived = {
        'B': (col['V'] + col['B-V'], [col['e_V'], col['e_B-V']]),
        'U': (col['V'] + col['B-V'] + col['U-B'], [col['e_V'], col['e_B-V'], col['e_U-B']]),
        'R': (col['V'] - col['V-R'], [col['e_V'], col['e_V-R']]),
        'I': (col['V'] - col['V-I'], [col['e_V'], col['e_V-I']])
    }
    for band, (mag, errs) in derived.items():
        err = np.ma.sqrt(sum(e**2 for e in errs))
        tbl[band] = MaskedColumn(mag.filled(0.), mask=np.ma.getmaskarray(mag), unit='mag')
        tbl[f'e_{band}'] = MaskedColumn(err.filled(0.), mask=np.ma.getmaskarray(err), unit='mag')


def build_catalog(ofile=None, match_radius=1.):
    r"""
    Construct the Landolt catalog distributed with pipnick from VizieR.

    This requires network access, and is only needed to regenerate the
    distributed catalog.  The VizieR catalogs ``J/AJ/137/4186`` (Landolt 2009)
    and ``J/AJ/146/131`` (Landolt 2013) are retrieved in full with
    :class:`astroquery.vizier.VizierClass`, standardized, and merged.  A star
    in both catalogs (within ``match_radius``) keeps its Landolt (2013) values.

    The merged catalog has the columns:

        - ``name``: the star name;
        - ``ra``, ``dec``: the J2000 coordinates (deg) at epoch 2000.0
          (``ref_epoch``);
        - ``pmra``, ``pmdec``: the proper motions (mas/yr; ``pmra`` includes
          the :math:`\cos\delta` factor);
        - ``V``, ``B-V``, ``U-B``, ``V-R``, ``R-I``, ``V-I`` and their errors
          (``e_V``, etc.), as published;
        - ``U``, ``B``, ``R``, ``I`` and their errors, derived from V and the
          colors; see :func:`derive_magnitudes`;
        - ``nobs``, ``nnight``: the numbers of observations and nights;
        - ``source``: the catalog reference.

    Parameters
    ----------
    ofile : str, :class:`pathlib.Path`, optional
        The output file.  If None, the catalog distributed with pipnick is
        overwritten; see :func:`catalog_path`.
    match_radius : float, optional
        The radius (arcsec) used to identify stars in both catalogs.

    Returns
    -------
    :class:`astropy.table.Table`
        The merged catalog.
    """
    vizier = Vizier(row_limit=-1, columns=['**'])
    common = {key: key for key in ['B-V', 'U-B', 'V-R', 'R-I', 'V-I']}
    common.update({f'e_{key}': f'e_{key}' for key in common})
    l09 = _ingest(
        vizier.get_catalogs('J/AJ/137/4186')['J/AJ/137/4186/stars'],
        {
            'name': 'Name', 'V': 'Vmag', 'e_V': 'e_Vmag', 'nobs': 'Nobs', 'nnight': 'Nnig',
            **common
        },
        'Landolt (2009, AJ, 137, 4186)'
    )
    l13_names = {key: f'<{key}>' for key in ['B-V', 'U-B', 'V-R', 'R-I', 'V-I']}
    l13_names.update({f'e_{key}': f'e_<{key}>' for key in l13_names})
    l13 = _ingest(
        vizier.get_catalogs('J/AJ/146/131')['J/AJ/146/131/standards'],
        {
            'name': 'Name', 'V': '<Vmag>', 'e_V': 'e_<Vmag>', 'nobs': 'No', 'nnight': 'Nn',
            **l13_names
        },
        'Landolt (2013, AJ, 146, 131)'
    )

    # Remove the Landolt (2009) stars that are also in Landolt (2013)
    c09 = SkyCoord(l09['ra'], l09['dec'], unit='deg')
    c13 = SkyCoord(l13['ra'], l13['dec'], unit='deg')
    _, sep, _ = c09.match_to_catalog_sky(c13)
    duplicate = sep.arcsec < match_radius
    if np.any(duplicate):
        log.info(
            f'{np.sum(duplicate)} Landolt (2009) stars are also in Landolt (2013); keeping the '
            'Landolt (2013) values.'
        )
    tbl = vstack([l09[np.logical_not(duplicate)], l13])
    derive_magnitudes(tbl)
    # Round the values to their published precision, so that the conversion
    # from single precision does not add spurious digits
    for name in tbl.colnames:
        if tbl[name].dtype.kind == 'f' and name != 'ref_epoch':
            tbl[name] = np.round(tbl[name], 7 if name in ['ra', 'dec'] else 4)
    tbl.meta = {
        'description': 'Landolt UBVRI standard stars, merged from Landolt (2009) and (2013).',
        'references': [
            'Landolt (2009, AJ, 137, 4186); VizieR J/AJ/137/4186',
            'Landolt (2013, AJ, 146, 131); VizieR J/AJ/146/131'
        ],
        'acknowledgment': 'This research has made use of the VizieR catalogue access tool, '
                           'CDS, Strasbourg, France (Ochsenbein et al. 2000, A&AS, 143, 23).',
        'duplicates_removed': int(np.sum(duplicate))
    }
    _ofile = catalog_path() if ofile is None else Path(ofile)
    _ofile.parent.mkdir(parents=True, exist_ok=True)
    tbl.write(_ofile, format='ascii.ecsv', overwrite=True)
    log.info(f'Wrote {len(tbl)} Landolt standards to {_ofile}')
    return tbl


def load_landolt(path=None):
    """
    Read the Landolt catalog.

    Parameters
    ----------
    path : str, :class:`pathlib.Path`, optional
        The catalog file.  If None, the catalog distributed with pipnick is
        read; see :func:`catalog_path`.

    Returns
    -------
    :class:`astropy.table.Table`
        The catalog; see :func:`build_catalog`.
    """
    return Table.read(catalog_path() if path is None else path, format='ascii.ecsv')


def catalog_color(catalog, color):
    """
    Return a color of the catalog stars and its error.

    Parameters
    ----------
    catalog : :class:`astropy.table.Table`
        The catalog; see :func:`build_catalog`.
    color : str
        The color, given by two bands separated by a hyphen, e.g. ``'B-V'``.
        The published colors (``B-V``, ``U-B``, ``V-R``, ``R-I``, ``V-I``) are
        used directly; any other color is the difference of the derived
        magnitudes.

    Returns
    -------
    value, error : :class:`numpy.ma.MaskedArray`
        The color and its error, masked where unavailable.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if the color is not a difference of two UBVRI bands.
    """
    bands = parse_color(color)
    if color in catalog.colnames:
        return (
            np.ma.asarray(catalog[color], dtype=float),
            np.ma.asarray(catalog[f'e_{color}'], dtype=float)
        )
    m1 = np.ma.asarray(catalog[bands[0]], dtype=float)
    m2 = np.ma.asarray(catalog[bands[1]], dtype=float)
    e1 = np.ma.asarray(catalog[f'e_{bands[0]}'], dtype=float)
    e2 = np.ma.asarray(catalog[f'e_{bands[1]}'], dtype=float)
    return m1 - m2, np.ma.sqrt(e1**2 + e2**2)


def parse_color(color):
    """
    Return the two bands of a color, e.g., ``('B', 'V')`` for ``'B-V'``.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if the color is not a difference of two different UBVRI bands.
    """
    bands = color.split('-')
    if (
        len(bands) != 2 or any(b not in 'UBVRI' or len(b) != 1 for b in bands)
        or bands[0] == bands[1]
    ):
        raise PipnickParameterError(
            f'Invalid color {color}; must be two different UBVRI bands, e.g. B-V.'
        )
    return bands[0], bands[1]


def match_landolt(catalog, ra, dec, epoch, radius=2.):
    """
    Match sources to the Landolt standards.

    The catalog positions are propagated to the epoch of the observation
    using their proper motions; see
    :func:`~pipnick.astrometry.services.propagate_proper_motion`.  Each
    standard is matched to its nearest source within ``radius``, and a source
    matched to more than one standard is kept only for the nearest; see
    :func:`~pipnick.utils.crossmatch.match_sky`.

    Parameters
    ----------
    catalog : :class:`astropy.table.Table`
        The Landolt catalog; see :func:`build_catalog`.
    ra, dec : :class:`numpy.ndarray`
        The ICRS coordinates of the sources (deg).
    epoch : float
        The epoch of the observation (Julian year).
    radius : float, optional
        The match radius (arcsec).

    Returns
    -------
    isrc, icat : :class:`numpy.ndarray`
        The indices of the matched sources and standards.
    sep : :class:`numpy.ndarray`
        The separations of the matches (arcsec).
    """
    if len(ra) == 0 or len(catalog) == 0:
        return np.array([], dtype=int), np.array([], dtype=int), np.array([])
    stars = propagate_proper_motion(catalog, epoch)
    sources = SkyCoord(np.asarray(ra, dtype=float), np.asarray(dec, dtype=float), unit='deg')
    icat, isrc, sep = match_sky(stars, sources, radius, unique='second')
    return isrc, icat, sep
