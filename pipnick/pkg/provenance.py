"""
Hashes used to record the provenance of the pipnick products.

Each step records hashes of its inputs and settings with its outputs (in
header keywords or table metadata), and compares them on later runs to
decide whether an output is up to date.  Every hash is the first 16
hexadecimal digits of a SHA-1 digest (see :func:`digest`).

.. warning::

    The hashes are stored in existing products, so the value each function
    returns for a given input must never change; otherwise, every product
    would be considered out of date and rebuilt.  The tests check the values
    against fixed references.
"""

import hashlib
import json

from pathlib import Path


def digest(*parts):
    """
    Return the hash of a sequence of strings or bytes.

    Parameters
    ----------
    *parts : str, bytes
        The content to hash, in order.  Strings are encoded as UTF-8.  The
        hash of several parts is the hash of their concatenation.

    Returns
    -------
    str
        The first 16 hexadecimal digits of the SHA-1 digest.
    """
    sha = hashlib.sha1()
    for part in parts:
        sha.update(part.encode() if isinstance(part, str) else part)
    return sha.hexdigest()[:16]


def json_hash(content):
    """
    Return the hash of the JSON representation of an object.

    Parameters
    ----------
    content : object
        Any object that can be serialized by :func:`json.dumps`, e.g., a
        dictionary of settings.  The keys of dictionaries are sorted, so the
        hash does not depend on their order.

    Returns
    -------
    str
        The hash; see :func:`digest`.
    """
    return digest(json.dumps(content, sort_keys=True))


def file_hash(path):
    """
    Return the hash of the content of a file.

    Parameters
    ----------
    path : str, :class:`pathlib.Path`
        The file.

    Returns
    -------
    str
        The hash; see :func:`digest`.
    """
    return digest(Path(path).read_bytes())


def names_hash(files):
    """
    Return the hash of the names of a set of files.

    Parameters
    ----------
    files : list
        The files (str or :class:`pathlib.Path`).  Only their names, without
        their directories, are used, and their order does not matter.

    Returns
    -------
    str
        The hash of the sorted names, separated by newlines; see
        :func:`digest`.
    """
    return digest('\n'.join(sorted(Path(f).name for f in files)))


def values_hash(values):
    """
    Return the hash of a set of named values.

    Parameters
    ----------
    values : dict
        The values, which are converted to strings with :func:`str`.  The
        order of the keys does not matter.

    Returns
    -------
    str
        The hash of the sorted ``key=value`` pairs, separated by newlines; see
        :func:`digest`.
    """
    return digest('\n'.join(f'{k}={values[k]}' for k in sorted(values.keys())))


def wcs_hash(wcs):
    """
    Return the hash of a WCS, used to determine if an astrometric solution has
    changed.

    Parameters
    ----------
    wcs : :class:`astropy.wcs.WCS`
        The WCS.

    Returns
    -------
    str
        The hash of its FITS header (with SIP distortion, if any); see
        :func:`digest`.
    """
    return digest(wcs.to_header(relax=True).tostring())


def image_hash(header):
    """
    Return the hash of the content of a reduced image.

    The hash combines the provenance keywords of the image (``PNRAWHSH`` and
    ``PNTABH``; see :func:`~pipnick.pipelines.reduction.reduce_all`), so it
    changes when the image is reduced from different frames or table values.

    Parameters
    ----------
    header : :class:`astropy.io.fits.Header`
        The header of the reduced image.

    Returns
    -------
    str
        The hash, ``<PNRAWHSH>-<PNTABH>``, with ``none`` for missing keywords.
    """
    return f'{header.get("PNRAWHSH", "none")}-{header.get("PNTABH", "none")}'


def solution_hash(header):
    """
    Return the hash of the astrometric solution of a reduced image.

    Only a solution found by pipnick (with the ``PNASTMTH`` keyword) counts;
    see :func:`~pipnick.pipelines.astrometry.apply_wcs`.

    Parameters
    ----------
    header : :class:`astropy.io.fits.Header`
        The header of the reduced image.

    Returns
    -------
    str
        The hash of the solution (``PNASTWH``), or ``'none'`` if the image has
        no solution.
    """
    return header.get('PNASTWH', 'none') if 'PNASTMTH' in header else 'none'
