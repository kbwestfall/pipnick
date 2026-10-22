"""
Basic image processing of Nickel data: overscan subtraction, trimming, bias
subtraction, field-flattening, and cosmic-ray masking, and the construction of
the combined calibration frames.

Frames are processed in the combination groups defined by the ``comb_id``
column of a per-setup metadata file; see :mod:`pipnick.metadata` and
:func:`reduce_all`.
"""

from pathlib import Path
import hashlib
import re
import warnings

from IPython import embed

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.wcs.wcs import FITSFixedWarning
from astropy import units
import ccdproc

from pipnick import cameras
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickPathError
from pipnick.metadata import build_metadata, frame_path, setup_tag, update_metadata_meta
from pipnick import log


def reduce_all(table, rdxdir=None, overwrite=False):
    """
    Perform reduction of raw astronomical data frames (overscan subtraction,
    bias subtraction, flat division, cosmic ray masking).

    The processing steps use :func:`ccdproc.subtract_overscan`,
    :func:`ccdproc.trim_image`, :func:`ccdproc.subtract_bias`,
    :func:`ccdproc.flat_correct`, and :func:`ccdproc.cosmicray_lacosmic`.

    The frames are reduced in the combination groups set by the ``comb_id``
    column of the metadata table; see :func:`check_comb_groups` for the rules
    the groups must follow.  Each group produces one output file in ``rdxdir``:

        - ``Bias.fits``: the combined bias.
        - ``Dark_<exptime>s.fits``: the combined, bias-subtracted dark for
          each exposure time (rounded to 0.1 s), e.g., ``Dark_120s.fits``.
          The darks are not yet applied to the science frames.
        - ``SkyFlat_<filter>.fits`` and ``DomeFlat_<filter>.fits``: the
          combined, bias-subtracted flats for each filter.
        - ``<object>_<filter>_C<comb_id>_<date>.fits``: the processed science
          frames, e.g., ``M11_V_C002_2026-08-19.fits``, where ``<object>`` is
          the most common object name in the group, ``<filter>`` is the
          filter (with a prime written as ``p``, e.g., ``ip`` for ``i'``),
          ``<comb_id>`` is the zero-padded group number, and ``<date>`` is the
          UT date of the group's first frame.  Science frames are
          bias-subtracted, flat-fielded with the sky flat for their filter
          (or the dome flat, if there is no sky flat), and cleaned of cosmic
          rays.  The frames in a group are combined without registration,
          after scaling them to the median exposure time of the group.

    The metadata table, not the raw headers, is the authority for the object,
    filter, exposure time, and coordinates, so that edits to the table (e.g.,
    a corrected object name) propagate to the outputs.  These values are
    written to the output headers, and any raw header value that differs is
    kept in an ``ORIG*`` keyword (e.g., ``ORIGOBJ``); see
    :func:`_apply_table_values`.

    The header of each output records the setup, the combination group, and
    the raw frames and calibrations used to construct it (as HISTORY lines and
    the ``PN*`` keywords).  When the reduction is rerun, an existing output is
    only reprocessed if its raw frames, calibrations, or table values (above)
    have changed, so that
    new frames can be added to the metadata table (see ``pipnick_metadata
    --split --append``) and the reduction simply rerun.  Outputs for the same
    setup that are no longer produced, e.g., because their group was merged
    into another, are reported but not deleted.

    A reduction directory may only hold the products from a single setup.

    Parameters
    ----------
    table : str, :class:`pathlib.Path`
        A per-setup metadata file, written by ``pipnick_metadata --split``;
        see :func:`~pipnick.metadata.build_metadata`.
    rdxdir : str, :class:`pathlib.Path`, optional
        Directory for the reduced data.  If None, the directory saved in the
        metadata file is used or, if there is none, the parent directory of the
        raw data directory.  I.e., if the raw directory is
        ``/User/janedoe/Nickel/raw``, ``rdxdir`` will be set to
        ``/User/janedoe/Nickel/``.  The directory is saved to the metadata
        file.
    overwrite : bool, optional
        Reprocess all outputs, even if they are up to date.  The metadata
        ``table`` is never rebuilt, to preserve any user edits; use
        ``pipnick_metadata --overwrite`` to rebuild it.

    Returns
    -------
    list
        Paths to the reduced science images, whether reprocessed or already up
        to date.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickPathError`
        Raised if the metadata file does not exist.
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the table includes more than one setup, has no frames, has
        invalid combination groups, or if the reduction directory includes
        products from a different setup.
    """
    _table = Path(table).absolute()
    if not _table.is_file():
        raise PipnickPathError(f'{_table} does not exist.')
    # NOTE: overwrite applies to the reduction products, not the metadata table;
    # rebuilding an existing table would discard any user edits to it.
    metafile, metadata = build_metadata(filename=_table, rdxdir=rdxdir)
    check_single_setup(metadata, metafile)
    if len(metadata) == 0:
        raise PipnickDataModelError(f'{metafile} has no frames to reduce.')
    check_comb_groups(metadata, metafile)
    tag = setup_tag(metadata.meta['setups'][metadata['setup'][0]])
    camera = cameras.get_camera(metadata.meta['camera'])

    # Get the output directory
    if 'rdxdir' in metadata.meta:
        _rdxdir = Path(metadata.meta['rdxdir'])
    else:
        _rdxdir = Path(metadata.meta['rawdir']).absolute().parent
        metadata.meta['rdxdir'] = str(_rdxdir)
        # Update the file to include the reduction directory
        update_metadata_meta(metafile, rdxdir=str(_rdxdir))
    if not _rdxdir.is_dir():
        _rdxdir.mkdir(parents=True)

    log.info(f'Reducing data:')
    log.info(f'    Metadata file: {metafile}')
    log.info(f'    Setup: {tag}')
    log.info(f'    Reduction directory: {_rdxdir}')

    # Set the output files and check that the directory does not include
    # products from another setup
    products = _plan_products(metadata, _rdxdir)
    existing = _find_products(_rdxdir)
    other = [s for s in existing.values() if s != tag]
    if len(other) > 0:
        raise PipnickDataModelError(
            f'{_rdxdir} includes {len(other)} products from a different setup than {tag} '
            f'({", ".join(sorted(set(other)))}).  A reduction directory can only hold the '
            'products from one setup; use a different directory.'
        )

    raw = frame_path(metadata, path_key='rawdir')

    # Combine the biases
    bias = None
    for product in products:
        if product['frametype'] == 'bias':
            bias = _reduce_calibration(product, metadata, raw, tag, overwrite, camera)
    if bias is None:
        log.warning('No bias frames available.  Attempting to continue anyway...')

    # Combine the darks and flats, which depend on the bias
    flats = {}
    for product in products:
        if product['frametype'] == 'dark':
            _reduce_calibration(product, metadata, raw, tag, overwrite, camera, bias=bias)
        elif product['frametype'] in ['domeflat', 'skyflat']:
            flats[(product['frametype'], product['filter'])] = _reduce_calibration(
                product, metadata, raw, tag, overwrite, camera, bias=bias
            )

    # Reduce the science frames, which depend on the bias and flat
    science = []
    for product in products:
        if product['frametype'] != 'science':
            continue
        flat = _select_flat(flats, product['filter'])
        _reduce_science(product, metadata, raw, tag, overwrite, camera, bias=bias, flat=flat)
        science.append(product['file'])
    if len(science) == 0:
        log.warning('No science frames available.')

    # Report any outputs for this setup that are no longer produced
    planned = [p['file'] for p in products]
    stale = sorted(
        f.name for f, s in _find_products(_rdxdir).items() if s == tag and f not in planned
    )
    if len(stale) > 0:
        log.warning(
            f'These files in {_rdxdir} are from setup {tag} but were not produced by this '
            f'reduction and may be out of date: {", ".join(stale)}.  They have not been deleted.'
        )

    return science


def check_single_setup(metadata, metafile):
    """
    Check that a metadata table includes frames from only one setup.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table, excluding commented-out rows; see
        :func:`~pipnick.metadata.read_metadata`.
    metafile : str, :class:`pathlib.Path`
        The file the table was read from, used in error messages.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the table has no ``setup`` column, or includes frames from
        more than one setup.
    """
    if 'setup' not in metadata.colnames:
        raise PipnickDataModelError(
            f'{metafile} has no setup column; it was likely built by an older version of '
            'pipnick.  Rebuild it using pipnick_metadata.'
        )
    setups, counts = np.unique(metadata['setup'], return_counts=True)
    if len(setups) > 1:
        definitions = metadata.meta.get('setups', {})
        summary = '\n'.join(
            f'    {s} ({n} frames): {definitions.get(s)}' for s, n in zip(setups, counts)
        )
        raise PipnickDataModelError(
            f'{metafile} includes frames from more than one setup:\n{summary}\nData from only '
            'one setup can be reduced at a time.  Comment out the rows for the other setups, '
            'or use pipnick_metadata --split to write one file per setup.'
        )


def check_comb_groups(metadata, metafile):
    """
    Check that the combination groups of a metadata table can be reduced.

    The groups are set by the ``comb_id`` column.  So that the calibrations
    used for each frame are unambiguous, the groups must follow these rules:

        - all frames in a group have the same frame type;
        - there is only one group of biases;
        - the frames in a group of darks have the same exposure time (to
          within 0.1 s), and there is only one group of darks for each
          exposure time;
        - the frames in a group of dome or sky flats have the same filter, and
          there is only one group of each flat type for each filter;
        - the frames in a group of science frames have the same filter.

    Science frames in a group with different exposure times are allowed, but a
    warning is issued.  Frames with other types (e.g., focus or pointing
    frames) are not reduced, and are only noted.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table, excluding commented-out rows; see
        :func:`~pipnick.metadata.read_metadata`.
    metafile : str, :class:`pathlib.Path`
        The file the table was read from, used in error messages.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the table has no ``comb_id`` column, or if any of the rules
        above are broken.  All broken rules are listed.
    """
    if 'comb_id' not in metadata.colnames:
        raise PipnickDataModelError(
            f'{metafile} has no comb_id column.  The reduction requires a per-setup metadata '
            'file, which defines the frames to combine; use pipnick_metadata --split to write '
            'one.'
        )
    violations = []
    bias_groups = []
    darks = {}
    flats = {}
    ignored = []
    for comb_id, indx in _group_indices(metadata).items():
        files = ', '.join(metadata['file'][indx])
        frametypes = sorted(set(metadata['frametype'][indx]))
        if len(frametypes) > 1:
            violations.append(
                f'Group {comb_id} mixes frame types ({", ".join(frametypes)}): {files}'
            )
            continue
        frametype = frametypes[0]
        filters = sorted(set(_filter(metadata, i) for i in indx))
        exptimes = set(_exptime(metadata, i) for i in indx)
        exptimes = set(None if e is None else round(e, 1) for e in exptimes)
        if frametype == 'bias':
            bias_groups.append(comb_id)
        elif frametype == 'dark':
            if len(exptimes) > 1:
                violations.append(
                    f'Dark group {comb_id} mixes exposure times '
                    f'({", ".join(sorted(str(e) for e in exptimes))}): {files}'
                )
            elif None in exptimes:
                violations.append(f'Dark group {comb_id} has no exposure time: {files}')
            else:
                darks.setdefault(exptimes.pop(), []).append(comb_id)
        elif frametype in ['domeflat', 'skyflat']:
            if len(filters) > 1:
                violations.append(
                    f'{frametype} group {comb_id} mixes filters ({", ".join(filters)}): {files}'
                )
            else:
                flats.setdefault((frametype, filters[0]), []).append(comb_id)
        elif frametype == 'science':
            if len(filters) > 1:
                violations.append(
                    f'Science group {comb_id} mixes filters ({", ".join(filters)}): {files}'
                )
            if len(exptimes) > 1:
                log.warning(
                    f'Science group {comb_id} mixes exposure times; the frames will be scaled '
                    f'to the median exposure time before they are combined: {files}'
                )
        else:
            ignored.extend(metadata['file'][indx])

    if len(bias_groups) > 1:
        violations.append(
            f'There is more than one bias group ({", ".join(str(c) for c in bias_groups)})'
        )
    for exptime, ids in darks.items():
        if len(ids) > 1:
            violations.append(
                f'There is more than one group of {exptime:g} s darks '
                f'({", ".join(str(c) for c in ids)})'
            )
    for (frametype, filt), ids in flats.items():
        if len(ids) > 1:
            violations.append(
                f'There is more than one group of {filt} {frametype} frames '
                f'({", ".join(str(c) for c in ids)})'
            )
    if len(ignored) > 0:
        log.info(f'These frames have types that are not reduced: {", ".join(ignored)}')
    if len(violations) > 0:
        raise PipnickDataModelError(
            f'{metafile} has invalid combination groups (comb_id):\n'
            + '\n'.join(f'    {v}' for v in violations)
        )


def _masked_to_none(value):
    """
    Return None if a table value is masked, or the value otherwise.
    """
    return None if value is np.ma.masked else value


def _filter(metadata, i):
    """
    Return the filter of a frame in the metadata table, or an empty string if
    it is unknown.
    """
    value = _masked_to_none(metadata['filter'][i])
    return '' if value is None else str(value)


def _exptime(metadata, i):
    """
    Return the exposure time of a frame in the metadata table, or None if it is
    unknown.
    """
    value = _masked_to_none(metadata['exptime'][i])
    return None if value is None else float(value)


def _group_indices(metadata):
    """
    Find the rows of the metadata table in each combination group.

    Returns
    -------
    dict
        The row indices of each group, keyed by the ``comb_id`` and ordered by
        it.  The indices are sorted by the MJD of the frame (frames without an
        MJD are last), and then by the row order.
    """
    mjd = np.ma.filled(metadata['mjd'], np.inf)
    srt = np.lexsort((np.arange(len(metadata)), mjd))
    comb_id = np.asarray(metadata['comb_id'])
    return {int(c): srt[comb_id[srt] == c] for c in np.unique(comb_id)}


def _safe_name(value):
    """
    Convert a string for use in a file name.

    Spaces are replaced by underscores, and any character other than letters,
    digits, underscores, periods, plus signs, and minus signs is removed.  For
    example, ``NGC 7331 I`` becomes ``NGC_7331_I``.  See :func:`_filter_name`
    for filters.
    """
    return re.sub(r'[^A-Za-z0-9_.+-]', '', str(value).strip().replace(' ', '_'))


def _filter_name(filt):
    """
    Convert a filter name for use in a file name.

    A prime is replaced by ``p``, so that primed filters are distinct from the
    unprimed filters with the same letter; e.g., ``i'`` becomes ``ip``, while
    ``i`` is unchanged.  The result is then converted by :func:`_safe_name`, or
    is ``none`` if the filter is unknown.
    """
    name = _safe_name(str(filt).replace("'", 'p'))
    return 'none' if len(name) == 0 else name


def _dark_name(exptime):
    """
    Return the name of the combined dark for an exposure time, which is rounded
    to 0.1 s; e.g., ``Dark_120s.fits`` or ``Dark_4.3s.fits``.
    """
    return f'Dark_{round(exptime, 1):g}s.fits'


def _flat_name(frametype, filt):
    """
    Return the name of the combined flat of a given type (``'skyflat'`` or
    ``'domeflat'``) and filter; e.g., ``SkyFlat_V.fits`` or
    ``SkyFlat_ip.fits`` for the ``i'`` filter (see :func:`_filter_name`).
    """
    title = 'SkyFlat' if frametype == 'skyflat' else 'DomeFlat'
    return f'{title}_{_filter_name(filt)}.fits'


def _group_object(metadata, indx):
    """
    Return the most common object name in a combination group.

    Ties go to the earliest frame.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table.
    indx : :class:`numpy.ndarray`
        The rows in the group, sorted by MJD; see :func:`_group_indices`.

    Returns
    -------
    str, None
        The object name, as given in the table, or None if it is unknown.
    """
    objects = [_masked_to_none(metadata['object'][i]) for i in indx]
    # max returns the first of the most common objects, i.e., the earliest
    return max(objects, key=objects.count)


def _table_values(metadata, indx, frametype):
    """
    Return the metadata table values to write to the header of an output.

    The table is the authority for these values, so that edits to it (e.g., a
    corrected object name) propagate to the reduced images.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table.
    indx : :class:`numpy.ndarray`
        The rows in the group, sorted by MJD; see :func:`_group_indices`.
    frametype : str
        The frame type of the group.

    Returns
    -------
    dict
        The values keyed by the quantity (``'filter'`` and ``'exptime'`` for
        all frame types, and ``'object'``, ``'ra'``, and ``'dec'`` for science
        frames).  Unknown values are None.  The object is the most common in
        the group (see :func:`_group_object`), the exposure time is the
        median, and the other values are from the earliest frame.
    """
    exptimes = [_exptime(metadata, i) for i in indx]
    values = {
        'filter': _masked_to_none(metadata['filter'][indx[0]]),
        'exptime': None if None in exptimes else float(np.median(exptimes))
    }
    if frametype == 'science':
        values['object'] = _group_object(metadata, indx)
        values['ra'] = _masked_to_none(metadata['ra'][indx[0]])
        values['dec'] = _masked_to_none(metadata['dec'][indx[0]])
    return values


def _table_hash(values):
    """
    Return a hash of the metadata table values written to the header of an
    output (see :func:`_table_values`), used to determine if they have changed.

    The hash is the first 16 hexadecimal digits of the SHA-1 hash of the
    sorted ``quantity=value`` pairs, separated by newlines.
    """
    text = '\n'.join(f'{k}={values[k]}' for k in sorted(values.keys()))
    return hashlib.sha1(text.encode()).hexdigest()[:16]


def _apply_table_values(ccd, values, camera):
    """
    Write metadata table values to the header of an output.

    Each value is written to the camera keyword for its quantity found in the
    header (e.g., ``FILTNAM`` for the filter of the original Nickel camera),
    and to the standard keyword (``OBJECT``, ``FILTER``, ``EXPTIME``, ``RA``,
    or ``DEC``) if that is different.  If the header value differs from the
    table value, the header value is saved to ``ORIGOBJ``, ``ORIGFILT``,
    ``ORIGEXP``, ``ORIGRA``, or ``ORIGDEC``.  Unknown (None) table values are
    not written.

    Parameters
    ----------
    ccd : :class:`astropy.nddata.CCDData`
        The output image, with a :class:`astropy.io.fits.Header` as its
        ``meta``.  Modified in place.
    values : dict
        The table values; see :func:`_table_values`.
    camera : :class:`~pipnick.cameras.cameras.NickelCamera`
        The camera used to obtain the data.
    """
    keys = {
        'object': ('OBJECT', 'ORIGOBJ'), 'filter': ('FILTER', 'ORIGFILT'),
        'exptime': ('EXPTIME', 'ORIGEXP'), 'ra': ('RA', 'ORIGRA'), 'dec': ('DEC', 'ORIGDEC')
    }
    for quantity, value in values.items():
        if value is None:
            continue
        standard, orig = keys[quantity]
        key = camera.header_keyword(ccd.meta, quantity)
        original = None if key is None else ccd.meta[key]
        if original is not None:
            if isinstance(value, float):
                changed = not np.isclose(float(original), value)
            else:
                changed = str(original).strip() != str(value).strip()
            if changed:
                ccd.meta[orig] = (original, f'{key} in the raw header')
        for k in ([standard] if key is None or key == standard else [key, standard]):
            ccd.meta[k] = value


def _science_name(metadata, indx, comb_id):
    """
    Return the name of the output science image for a combination group.

    The name is ``<object>_<filter>_C<comb_id>_<date>.fits``, where
    ``<object>`` is the most common object name in the group (ties go to the
    earliest frame), ``<filter>`` is the filter of the earliest frame (all
    frames in a group have the same filter; see :func:`check_comb_groups` and
    :func:`_filter_name`), ``<comb_id>`` is zero-padded to three digits, and
    ``<date>`` is the UT date of the earliest frame.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table.
    indx : :class:`numpy.ndarray`
        The rows in the group, sorted by MJD; see :func:`_group_indices`.
    comb_id : int
        The combination group.

    Returns
    -------
    str
        The file name.
    """
    obj = _group_object(metadata, indx)
    obj = '' if obj is None else _safe_name(obj)
    date = _masked_to_none(metadata['date'][indx[0]])
    date = 'nodate' if date is None else str(date).split('T')[0]
    filt = _filter_name(_filter(metadata, indx[0]))
    return f'{"unknown" if len(obj) == 0 else obj}_{filt}_C{comb_id:03d}_{date}.fits'


def _plan_products(metadata, rdxdir):
    """
    Set the output file for each combination group to be reduced.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table, with valid combination groups; see
        :func:`check_comb_groups`.
    rdxdir : :class:`pathlib.Path`
        The reduction directory.

    Returns
    -------
    list
        A dictionary for each group to be reduced, ordered by ``comb_id``, with
        the frame type, the ``comb_id``, the row indices of its frames (see
        :func:`_group_indices`), its filter and exposure time (from its first
        frame), and the path to its output file.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if two groups would have the same output file, e.g., for flats
        with filter names that only differ in characters that cannot be used
        in a file name.
    """
    products = []
    for comb_id, indx in _group_indices(metadata).items():
        frametype = metadata['frametype'][indx[0]]
        product = {
            'frametype': frametype, 'comb_id': comb_id, 'indx': indx,
            'filter': _filter(metadata, indx[0]), 'exptime': _exptime(metadata, indx[0])
        }
        if frametype == 'bias':
            name = 'Bias.fits'
        elif frametype == 'dark':
            name = _dark_name(product['exptime'])
        elif frametype in ['domeflat', 'skyflat']:
            name = _flat_name(frametype, product['filter'])
        elif frametype == 'science':
            name = _science_name(metadata, indx, comb_id)
        else:
            continue
        product['file'] = rdxdir / name
        products.append(product)
    names = [p['file'].name for p in products]
    duplicates = sorted(set(n for n in names if names.count(n) > 1))
    if len(duplicates) > 0:
        raise PipnickDataModelError(
            f'More than one combination group would be written to: {", ".join(duplicates)}.'
        )
    return products


def _find_products(rdxdir):
    """
    Find the pipnick products in a reduction directory.

    Parameters
    ----------
    rdxdir : :class:`pathlib.Path`
        The reduction directory.

    Returns
    -------
    dict
        The setup tag (the ``PNSETUP`` header keyword) of each FITS file in the
        directory that has one, keyed by its path.
    """
    products = {}
    for f in sorted(rdxdir.glob('*.fits')):
        try:
            hdr = fits.getheader(f)
        except OSError:
            continue
        if 'PNSETUP' in hdr:
            products[f] = hdr['PNSETUP']
    return products


def _raw_hash(files):
    """
    Return a hash of the names of a set of raw files, used to determine if the
    frames used to construct an output have changed.

    The hash is the first 16 hexadecimal digits of the SHA-1 hash of the
    sorted file names (without their paths), separated by newlines.
    """
    names = '\n'.join(sorted(Path(f).name for f in files))
    return hashlib.sha1(names.encode()).hexdigest()[:16]


def _needs_update(ofile, tag, raw_hash, table_hash, calib_hashes, overwrite):
    """
    Determine if an output file must be (re)constructed.

    Parameters
    ----------
    ofile : :class:`pathlib.Path`
        The output file.
    tag : str
        The setup tag; see :func:`~pipnick.metadata.setup_tag`.
    raw_hash : str
        The hash of the raw frames used to construct the output; see
        :func:`_raw_hash`.
    table_hash : str
        The hash of the metadata table values written to the output header;
        see :func:`_table_hash`.  Compared to the ``PNTABH`` header keyword.
    calib_hashes : dict
        The hash of each calibration used to construct the output, keyed by
        the calibration type (e.g., ``'bias'``), or ``'none'`` if the
        calibration is not available.  Each is compared to the ``PN<TYPE>H``
        header keyword, e.g., ``PNBIASH``.
    overwrite : bool
        Reconstruct the output even if it is up to date.

    Returns
    -------
    update : bool
        True if the output must be (re)constructed.
    reason : str
        The reason.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the existing file is from a different setup.
    """
    if not ofile.is_file():
        return True, 'new output'
    hdr = fits.getheader(ofile)
    if 'PNSETUP' not in hdr:
        return True, 'the existing file has no provenance keywords'
    if hdr['PNSETUP'] != tag:
        raise PipnickDataModelError(
            f'{ofile} is from setup {hdr["PNSETUP"]}, not {tag}.  A reduction directory can '
            'only hold the products from one setup.'
        )
    if overwrite:
        return True, 'overwrite is set'
    if hdr.get('PNRAWHSH') != raw_hash:
        return True, 'the raw frames changed'
    if hdr.get('PNTABH') != table_hash:
        return True, 'the metadata table values changed'
    for name, value in calib_hashes.items():
        if hdr.get(f'PN{name.upper()}H') != value:
            return True, f'the {name} changed'
    return False, 'up to date'


def _calib_hashes(calibs):
    """
    Return the hash of each calibration, or ``'none'`` if it is not available.
    """
    return {name: 'none' if calib is None else calib['hash'] for name, calib in calibs.items()}


def _calibration_data(calib):
    """
    Return the image data for a calibration, reading it if necessary.

    Parameters
    ----------
    calib : dict, None
        The calibration, as returned by :func:`_reduce_calibration`.  If the
        data have not yet been read, they are read and saved to the
        dictionary.

    Returns
    -------
    :class:`astropy.nddata.CCDData`, None
        The calibration data, or None if ``calib`` is None.
    """
    if calib is None:
        return None
    if calib['ccd'] is None:
        calib['ccd'] = CCDData.read(calib['file'])
    return calib['ccd']


def _add_provenance(ccd, tag, comb_id, files, raw_hash, calibs):
    """
    Record the setup, raw frames, and calibrations used to construct an output.

    Parameters
    ----------
    ccd : :class:`astropy.nddata.CCDData`
        The output image.  Its header is modified in place.
    tag : str
        The setup tag; see :func:`~pipnick.metadata.setup_tag`.
    comb_id : int
        The combination group.
    files : list
        The raw files.
    raw_hash : str
        The hash of the raw files; see :func:`_raw_hash`.
    calibs : dict
        The calibrations used, keyed by their type (e.g., ``'bias'``), with
        None for any that are not available; see :func:`_reduce_calibration`.
    """
    if not isinstance(ccd.meta, fits.Header):
        ccd.meta = fits.Header(ccd.meta)
    ccd.meta['PNSETUP'] = (tag, 'pipnick setup')
    ccd.meta['PNCOMBID'] = (comb_id, 'pipnick combination group')
    ccd.meta['PNNRAW'] = (len(files), 'Number of raw frames')
    ccd.meta['PNRAWHSH'] = (raw_hash, 'Hash of the raw frame names')
    for name, calib in calibs.items():
        key = name.upper()
        ccd.meta[f'PN{key}'] = ('none' if calib is None else calib['file'].name, f'{name} used')
        ccd.meta[f'PN{key}H'] = (
            'none' if calib is None else calib['hash'], f'Hash of the {name} raw frame names'
        )
    for f in files:
        ccd.meta.add_history(f'pipnick raw frame: {Path(f).name}')
    for calib in calibs.values():
        if calib is not None:
            ccd.meta.add_history(f'pipnick calibration: {calib["file"].name}')


def _select_flat(flats, filt):
    """
    Select the flat for a filter: the sky flat if available, otherwise the dome
    flat.

    Parameters
    ----------
    flats : dict
        The flats, keyed by their type and filter; see
        :func:`_reduce_calibration`.
    filt : str
        The filter.

    Returns
    -------
    dict, None
        The flat, or None if there is no flat for the filter.
    """
    for frametype in ['skyflat', 'domeflat']:
        if (frametype, filt) in flats:
            return flats[(frametype, filt)]
    log.warning(f'No flat available for filter {filt}.  Science frames will not be flat-fielded.')
    return None


def _check_pointing(metadata, indx, camera):
    """
    Warn if the frames in a group have different pointings.

    The frames are combined without registration, so a warning is issued if
    the header coordinates of any frame are more than 3 binned pixels from
    those of the first frame.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table.
    indx : :class:`numpy.ndarray`
        The rows in the group.
    camera : :class:`~pipnick.cameras.cameras.NickelCamera`
        The camera used to obtain the data, which sets the plate scale.
    """
    if len(indx) < 2:
        return
    ra = [_masked_to_none(metadata['ra'][i]) for i in indx]
    dec = [_masked_to_none(metadata['dec'][i]) for i in indx]
    if None in ra or None in dec:
        log.warning('Cannot check the pointing of frames without coordinates.')
        return
    try:
        coo = SkyCoord(ra, dec, unit=(units.hourangle, units.deg))
    except ValueError:
        log.warning('Cannot parse the frame coordinates to check their pointing.')
        return
    sep = coo[0].separation(coo).to('arcsec').value
    nbin = max(int(metadata['xbin'][indx[0]]), int(metadata['ybin'][indx[0]]))
    tolerance = 3 * camera.plate_scale * nbin
    offset = sep > tolerance
    if np.any(offset):
        log.warning(
            f'These frames are more than {tolerance:.1f} arcsec from '
            f'{metadata["file"][indx[0]]}, but are combined without registration: '
            f'{", ".join(metadata["file"][indx[offset]])}'
        )


def init_ccddata(frame, bpm=None, saturation=None, gain=None, readnoise=None):
    """
    Read data from a provided file.
    
    This function also uses the header information to identify the camera used
    to acquire the data.  If the camera is successfully identified, the
    camera information is used to 

        - apply any known bad-pixel mask
        - mask saturated pixels
        - convert the units from DN to counts by multiplying by the camera gain

    The relevant keyword arguments can be used to supplement information not
    available for a given camera or when the camera cannot be identified.  Note
    that values provided via a keyword *always* take precedence.

    Parameters
    ----------
    frame : str, :class:`pathlib.Path`
        Path to the FITS file.

    Returns
    -------
    :class:`astropy.nddata.CCDData`
        Initialized and processed CCDData object.
    """
#    warnings.simplefilter('ignore', FITSFixedWarning)
    _frame = Path(frame).absolute()

    log.info(f'Initializing CCD imaging data in {_frame.name}.')
    ccd = CCDData.read(_frame, unit=units.adu)

    # Determine the camera used
    try:
        camera = cameras.identify_camera(ccd.header)
    except PipnickDataModelError:
        warnings.warn(f'Unable to determine camera used for {_frame.name}.')
        camera = None
        readout = None
    else:
        ccd.meta['camera'] = camera.name
        readout = camera.readout(ccd.header)
        if readout.placeholder:
            log.warning(
                f'The {camera.name} gain and read noise for the {readout.mode} readout '
                'are placeholders.'
            )

    # Mask bad pixels
    _bpm = camera.bpm(ccd.header) if bpm is None and camera is not None else bpm
    if _bpm is not None:
        ccd.mask = _bpm

    # Mask saturated pixels; saturation is in ADU
    _sat = readout.saturation if saturation is None and readout is not None else saturation
    if _sat is not None:
        if ccd.mask is None:
            ccd.mask = np.zeros(ccd.data.shape, dtype=bool)
        saturated = ccd.data >= _sat
        nsat = int(np.sum(saturated))
        if nsat > 0:
            # Bad pixels can read as saturated, so these are reported separately
            nbad = int(np.sum(saturated & ccd.mask))
            log.info(
                f'{_frame.name} has {nsat} saturated pixels (>= {_sat} ADU), {nbad} of which '
                'are known bad pixels.'
            )
        ccd.mask[saturated] = True

    # Get and save the RN (in electrons) and gain (in electrons/ADU)
    ccd.meta['PROCGAIN'] = readout.gain if gain is None and readout is not None else gain
    ccd.meta['PROCRN'] = (
        readout.readnoise if readnoise is None and readout is not None else readnoise
    )

    # Apply gain if available
    if ccd.meta['PROCGAIN'] is not None:
        ccd.data = ccd.data * ccd.meta['PROCGAIN']
        ccd.unit *= (units.electron / units.adu)

    # Done
    return ccd


def process_frame(raw_frame, bias=None, flat=None, flag_cosmics=True):
    """
    Perform basic image processing

    This includes overscan subtraction, trimming, bias subtraction, and
    field-flattening.  The mean overscan level (in ADU), the size of the
    trimmed image, and the number of pixels affected by cosmic rays are
    logged, and the overscan level is saved to the ``PNOSCAN`` header keyword.

    Parameters
    ----------
    raw_frame : str, :class:`pathlib.Path`
        Data file to process.
    flag_cosmics : bool, optional
        Identify and mask cosmic rays
    bias : :class:`astropy.nddata.CCDData`, optional
        Image to use for the bias subtraction
    flat : :class:`astropy.nddata.CCDData`, optional
        Image to use for the field-flattening

    Returns
    -------
    :class:`astropy.nddata.CCDData`
        Processed CCDData object with overscan subtracted and image trimmed.
    """
#    warnings.simplefilter("error", RuntimeWarning)
    ccd = init_ccddata(raw_frame)
    if 'camera' not in ccd.meta:
        log.warning('Unknown camera; cannot determine the overscan and data sections!')
    else:
        datasec, biassec = cameras.get_camera(ccd.meta['camera']).sections(ccd.header)
        if biassec is None:
            log.warning('Frame has no overscan region.')
        else:
            level = _masked_mean(ccdproc.trim_image(ccd, fits_section=biassec))
            if ccd.meta['PROCGAIN'] is not None:
                level /= ccd.meta['PROCGAIN']
            ccd.meta['PNOSCAN'] = (level, 'Mean overscan level (ADU)')
            log.info(f'Subtracting overscan; mean level is {level:.1f} ADU.')
            ccd = ccdproc.subtract_overscan(ccd, fits_section=biassec, overscan_axis=1)
        ccd = ccdproc.trim_image(ccd, fits_section=datasec)
        log.info(f'Trimmed image to {ccd.shape[1]} x {ccd.shape[0]} (columns x rows).')
    if bias is not None:
        log.info('Subtracting super bias')
        ccd = ccdproc.subtract_bias(ccd, bias)
    if flat is not None:
        log.info('Applying flat-field correction')
        ccd = ccdproc.flat_correct(ccd, flat)
    if flag_cosmics:
        nmasked = 0 if ccd.mask is None else int(np.sum(ccd.mask))
        ccd = ccdproc.cosmicray_lacosmic(
            ccd, gain_apply=False, gain=1.0, readnoise=ccd.meta['PROCRN'], verbose=False
        )
        ncosmic = 0 if ccd.mask is None else int(np.sum(ccd.mask)) - nmasked
        log.info(f'Removed cosmic rays from {Path(raw_frame).name}; {ncosmic} pixels affected.')
    return ccd


def _masked_mean(ccd):
    """
    Return the mean of the unmasked pixels in an image.
    """
    data = np.ma.MaskedArray(ccd.data, mask=ccd.mask)
    return float(np.ma.mean(data))


def process_frames(raw_frames, bias=None, flat=None, flag_cosmics=True):
    """
    Perform basic image processing of a list of frames; see
    :func:`process_frame`.

    Parameters
    ----------
    raw_frames : list
        List of file names to process.
    bias : :class:`astropy.nddata.CCDData`, optional
        Image to use for the bias subtraction
    flat : :class:`astropy.nddata.CCDData`, optional
        Image to use for the field-flattening
    flag_cosmics : bool, optional
        Identify and mask cosmic rays

    Returns
    -------
    list
        The processed :class:`astropy.nddata.CCDData` objects.
    """
    # TODO: For *many* images, this may be too memory intensive.
    return [
        process_frame(frame, bias=bias, flat=flat, flag_cosmics=flag_cosmics)
        for frame in raw_frames
    ]


def combine_frames(frames, scale=None):
    """
    Combine a list of processed frames using a sigma-clipped mean.

    Outliers are iteratively rejected at 5 sigma.  The header and WCS of the
    combined frame are those of the first frame, with ``NCOMBINE`` set to the
    number of frames.  A single frame is returned without combination (but
    with any scaling applied).

    Parameters
    ----------
    frames : list
        The :class:`astropy.nddata.CCDData` objects to combine.
    scale : array-like, optional
        Factors used to multiply each frame before combining them.  If None,
        the frames are not scaled.

    Returns
    -------
    :class:`astropy.nddata.CCDData`
        Combined CCDData object.
    """
    # NOTE: The arithmetic and the combination do not propagate the header
    meta = frames[0].meta.copy()
    wcs = frames[0].wcs
    _frames = frames if scale is None else [f.multiply(s) for s, f in zip(scale, frames)]

    if len(_frames) == 1:
        combined = _frames[0]
    else:
        log.info('Identifying outliers using sigma clipping')
        combiner = ccdproc.Combiner(_frames)

        old_n_masked = 0
        new_n_masked = 1
        while new_n_masked > old_n_masked:
            # TODO:
            #   - does ccdproc provide an iterative sigma clipper?
            #   - allow user to change the clipping parameters
            combiner.data_arr = sigma_clip(
                combiner.data_arr, sigma=5., maxiters=1, cenfunc=np.ma.mean, stdfunc=np.ma.std,
                axis=0, masked=True, copy=False
            )
            old_n_masked = new_n_masked
            new_n_masked = combiner.data_arr.mask.sum()

        log.info('Creating sigma-clipped average image.')
        combined = combiner.average_combine(scale_func=np.ma.mean)

    # TODO: Fill the masked regions?
    combined.meta = meta
    combined.meta['NCOMBINE'] = len(frames)
    combined.wcs = wcs
    return combined


def _report_bias_stats(frames, files):
    r"""
    Report the bias level of each bias frame and the read noise measured from
    the set of frames.

    The bias level is the mean of the unmasked pixels in the data section
    before overscan subtraction, in ADU, for comparison with the overscan
    level.  The read noise is measured from the difference of each consecutive
    pair of frames as :math:`\sigma_{\rm diff}/\sqrt{2}`, where
    :math:`\sigma_{\rm diff}` is the standard deviation of the difference after
    rejecting 5-sigma outliers (e.g., cosmic rays).  The frames include the
    nominal gain, so the read noise is in electrons (unless the gain is
    unknown).

    Parameters
    ----------
    frames : list
        The processed bias frames; see :func:`process_frames`.
    files : list
        The raw files.
    """
    for f, frame in zip(files, frames):
        gain = frame.meta['PROCGAIN']
        residual = _masked_mean(frame) / (1. if gain is None else gain)
        oscan = frame.meta.get('PNOSCAN', 0.)
        log.info(
            f'{Path(f).name}: mean bias level is {oscan + residual:.1f} ADU ({residual:.2f} ADU '
            'after overscan subtraction).'
        )
    if len(frames) < 2:
        log.info('At least two bias frames are needed to measure the read noise.')
        return

    readnoise = []
    for i in range(len(frames) - 1):
        mask = np.zeros(frames[i].shape, dtype=bool)
        for frame in frames[i:i+2]:
            if frame.mask is not None:
                mask |= frame.mask
        diff = np.ma.MaskedArray(frames[i+1].data - frames[i].data, mask=mask)
        _, _, std = sigma_clipped_stats(diff, sigma=5.)
        readnoise.append(std / np.sqrt(2.))
    unit = frames[0].unit.to_string()
    nominal = frames[0].meta['PROCRN']
    log.info(
        f'Read noise from {len(readnoise)} pairs of bias frames: median {np.median(readnoise):.2f}'
        f' {unit} (range {np.min(readnoise):.2f}-{np.max(readnoise):.2f}); nominal value is '
        f'{nominal} electron.'
    )


def _reduce_calibration(product, metadata, raw, tag, overwrite, camera, bias=None):
    """
    Construct a combined calibration frame, unless it is up to date.

    Biases are combined directly.  Darks and flats are bias-subtracted, and
    flats are scaled by their mean before they are combined.

    Parameters
    ----------
    product : dict
        The combination group; see :func:`_plan_products`.
    metadata : :class:`astropy.table.Table`
        The metadata table.
    raw : list
        The path to the raw file for each row in the metadata table.
    tag : str
        The setup tag; see :func:`~pipnick.metadata.setup_tag`.
    overwrite : bool
        Reconstruct the calibration even if it is up to date.
    camera : :class:`~pipnick.cameras.cameras.NickelCamera`
        The camera used to obtain the data.
    bias : dict, optional
        The bias, as returned by this function, used for darks and flats.

    Returns
    -------
    dict
        The calibration, with its output file (``file``), the hash of its raw
        frames (``hash``; see :func:`_raw_hash`), and its image data (``ccd``),
        which is None until needed if the calibration was up to date; see
        :func:`_calibration_data`.
    """
    frametype = product['frametype']
    files = [raw[i] for i in product['indx']]
    raw_hash = _raw_hash(files)
    calib = {'file': product['file'], 'hash': raw_hash, 'ccd': None}
    calibs = {} if frametype == 'bias' else {'bias': bias}
    values = _table_values(metadata, product['indx'], frametype)
    table_hash = _table_hash(values)
    update, reason = _needs_update(
        product['file'], tag, raw_hash, table_hash, _calib_hashes(calibs), overwrite
    )
    if not update:
        log.info(f'{product["file"].name} is up to date.')
        return calib

    log.info(f'Reprocessing {product["file"].name} ({reason}) from {len(files)} frames.')
    frames = process_frames(files, bias=_calibration_data(calibs.get('bias')), flag_cosmics=False)
    if frametype == 'bias':
        _report_bias_stats(frames, files)
    scale = None
    if frametype in ['domeflat', 'skyflat']:
        mean = np.array([float(f.mean().data) for f in frames])
        scale = mean[len(mean)//2] / mean
    ccd = combine_frames(frames, scale=scale)
    if frametype == 'bias':
        ccd.meta['OBJECT'] = 'Bias'
    elif frametype == 'dark':
        ccd.meta['OBJECT'] = 'Dark'
    else:
        ccd.meta['OBJECT'] = 'SkyFlat' if frametype == 'skyflat' else 'DomeFlat'
    _add_provenance(ccd, tag, product['comb_id'], files, raw_hash, calibs)
    _apply_table_values(ccd, values, camera)
    ccd.meta['PNTABH'] = (table_hash, 'Hash of the metadata table values')
    ccd.write(product['file'], overwrite=True)
    log.info(f'Saved {product["file"]}')
    calib['ccd'] = ccd
    return calib


def _reduce_science(product, metadata, raw, tag, overwrite, camera, bias=None, flat=None):
    """
    Process and combine a group of science frames, unless the output is up to
    date.

    The frames are bias-subtracted, flat-fielded, and cleaned of cosmic rays.
    If there is more than one frame, they are scaled to the median exposure
    time of the group and combined, without registration.  The object,
    filter, exposure time (the median), and coordinates in the header are set
    by the metadata table; see :func:`_apply_table_values`.

    Parameters
    ----------
    product : dict
        The combination group; see :func:`_plan_products`.
    metadata : :class:`astropy.table.Table`
        The metadata table.
    raw : list
        The path to the raw file for each row in the metadata table.
    tag : str
        The setup tag; see :func:`~pipnick.metadata.setup_tag`.
    overwrite : bool
        Reprocess the frames even if the output is up to date.
    camera : :class:`~pipnick.cameras.cameras.NickelCamera`
        The camera used to obtain the data.
    bias : dict, optional
        The bias; see :func:`_reduce_calibration`.
    flat : dict, optional
        The flat; see :func:`_reduce_calibration`.
    """
    indx = product['indx']
    files = [raw[i] for i in indx]
    raw_hash = _raw_hash(files)
    calibs = {'bias': bias, 'flat': flat}
    values = _table_values(metadata, indx, 'science')
    table_hash = _table_hash(values)
    update, reason = _needs_update(
        product['file'], tag, raw_hash, table_hash, _calib_hashes(calibs), overwrite
    )
    if not update:
        log.info(f'{product["file"].name} is up to date.')
        return

    log.info(f'Reprocessing {product["file"].name} ({reason}) from {len(files)} frames.')
    _check_pointing(metadata, indx, camera)
    frames = process_frames(
        files, bias=_calibration_data(bias), flat=_calibration_data(flat), flag_cosmics=True
    )
    scale = None
    if len(frames) > 1:
        if values['exptime'] is None:
            log.warning(
                'Some frames have no exposure time; the frames are combined without scaling.'
            )
        else:
            scale = values['exptime'] / np.array([_exptime(metadata, i) for i in indx])
    ccd = combine_frames(frames, scale=scale)
    _add_provenance(ccd, tag, product['comb_id'], files, raw_hash, calibs)
    _apply_table_values(ccd, values, camera)
    ccd.meta['PNTABH'] = (table_hash, 'Hash of the metadata table values')
    ccd.write(product['file'], overwrite=True)
    log.info(f'Saved {product["file"]}')
