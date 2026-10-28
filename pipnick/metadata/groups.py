"""
The groups of frames in a per-setup metadata file: the frames combined by
the reduction (``comb_id``), and the reduced science images combined into
mosaics (``mosaic`` and ``wcs_group``); see :mod:`pipnick.metadata`.
"""

import numpy as np
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.table import Column

from pipnick import cameras
from pipnick import log
from pipnick.metadata.names import filter_name, safe_name

mosaic_columns = ['mosaic', 'wcs_group']
"""
The columns of a per-setup file that define the mosaics, with string labels:
``mosaic`` identifies the reduced science images combined into one mosaic, and
``wcs_group`` the mosaics that share an output WCS (e.g., the mosaics of one
object in several filters).  An empty label means the frame is not
mosaicked.  The labels are global: the same label in the tables of different
setups (and, eventually, nights) refers to the same mosaic; see
:mod:`pipnick.pipelines.mosaic`.
"""


def _comb_key(row):
    """
    Return the key used to group a frame with others by default, or None if the
    frame is not combined with any other.

    Biases are combined, flats are combined by type and filter, and darks are
    combined if their exposure times agree to within 0.1 s; all within the same
    setup.  All other frames are processed independently.
    """
    frametype = row['frametype']
    if frametype == 'bias':
        return (row['setup'], frametype)
    if frametype in ['domeflat', 'skyflat']:
        filt = '' if row['filter'] is np.ma.masked else str(row['filter'])
        return (row['setup'], frametype, filt)
    if frametype == 'dark' and row['exptime'] is not np.ma.masked:
        return (row['setup'], frametype, round(float(row['exptime']), 1))
    return None



def assign_comb_ids(tbl, assign=None):
    """
    Assign the combination group of each frame in a metadata table.

    Frames with the same ``comb_id`` are combined.  By default, biases are
    combined, flats are combined by type and filter, darks are combined if
    their exposure times agree to within 0.1 s, and all other frames, including
    science frames, are processed independently; see :func:`_comb_key`.
    Groups never span setups.

    Parameters
    ----------
    tbl : :class:`astropy.table.Table`
        The metadata table, with the ``setup`` column.  The ``comb_id`` column
        is added after the ``frametype`` column if it does not exist.  Modified
        in place.
    assign : :class:`numpy.ndarray`, optional
        Boolean array selecting the frames to assign, e.g., newly appended
        frames.  The other frames keep their groups, and a new frame joins the
        existing group with the same default key, if there is exactly one.  If
        None, all frames are assigned.
    """
    _assign = np.ones(len(tbl), dtype=bool) if assign is None else assign
    ids = np.zeros(len(tbl), dtype=int)
    existing = {}
    if 'comb_id' in tbl.colnames:
        keep = np.logical_not(_assign)
        ids[keep] = np.asarray(tbl['comb_id'])[keep]
        for row, comb_id in zip(tbl[keep], ids[keep]):
            key = _comb_key(row)
            if key is not None:
                existing.setdefault(key, set()).add(int(comb_id))
    next_id = 1 if np.all(_assign) else int(np.max(ids[np.logical_not(_assign)])) + 1
    assigned = {}
    for i in np.where(_assign)[0]:
        key = _comb_key(tbl[i])
        if key is not None and key in assigned:
            ids[i] = assigned[key]
            continue
        if key is not None and key in existing and len(existing[key]) == 1:
            ids[i] = next(iter(existing[key]))
        else:
            ids[i] = next_id
            next_id += 1
        if key is not None:
            assigned[key] = ids[i]
    if 'comb_id' in tbl.colnames:
        tbl['comb_id'] = ids
    else:
        tbl.add_column(Column(ids, name='comb_id'), index=tbl.colnames.index('frametype') + 1)



def add_mosaic_columns(tbl):
    """
    Add the mosaic columns (see :data:`mosaic_columns`) to a metadata table,
    if they do not exist, and replace missing labels with empty strings.

    The columns are added after the ``comb_id`` column, if it exists, or at
    the end of the table.

    Parameters
    ----------
    tbl : :class:`astropy.table.Table`
        The metadata table.  Modified in place.
    """
    for name in mosaic_columns:
        if name in tbl.colnames:
            values = ['' if v is np.ma.masked else str(v) for v in tbl[name]]
            tbl[name] = Column(np.array(values + [''], dtype=str)[:len(values)], name=name)
            continue
        if name == 'wcs_group' and 'mosaic' in tbl.colnames:
            index = tbl.colnames.index('mosaic') + 1
        elif 'comb_id' in tbl.colnames:
            index = tbl.colnames.index('comb_id') + 1
        else:
            index = len(tbl.colnames)
        tbl.add_column(Column(np.full(len(tbl), '', dtype='U1'), name=name), index=index)



def _strip_filter(obj, filt):
    """
    Remove a trailing filter name from an object name, e.g. ``'M11'`` from
    ``'M11 V'``; observers often append the filter to the object name.

    The last word is removed if it is the filter of the frame, or a broad-band
    filter name (UBVRI or ugriz, with or without a prime).
    """
    words = str(obj).split()
    filters = {str(filt).strip()} | set('UBVRI') | {f'{b}{p}' for b in 'ugriz' for p in ['', "'"]}
    if len(words) > 1 and words[-1] in filters:
        words = words[:-1]
    return ' '.join(words)



def _friends_of_friends(coords, radius):
    """
    Group positions that are linked by separations smaller than ``radius``.

    Parameters
    ----------
    coords : :class:`astropy.coordinates.SkyCoord`
        The positions.
    radius : :class:`astropy.units.Quantity`
        The linking length.

    Returns
    -------
    :class:`numpy.ndarray`
        The group index of each position, numbered in order of first
        appearance.
    """
    n = len(coords)
    group = np.arange(n)
    for i in range(n):
        near = np.where(coords[i].separation(coords) < radius)[0]
        # Merge every group linked to position i
        linked = np.unique(group[near])
        group[np.isin(group, linked)] = linked[0]
    _, first, inverse = np.unique(group, return_index=True, return_inverse=True)
    order = np.argsort(np.argsort(first))
    return order[inverse]



def assign_mosaics(tbl, commented=None, radius=None):
    """
    Assign default mosaic labels to the science frames of a per-setup table.

    Only uncommented science frames without labels are assigned, so existing
    labels (e.g., user edits) are kept.  Combination groups (``comb_id``) are
    assigned as a whole, from the coordinates of their earliest frame.

    The frames are grouped by pointing: pointings within ``radius`` of each
    other, directly or through other pointings, form one group (a
    friends-of-friends grouping), which is given one WCS group label.  The
    object names are not used for the grouping, because they are often wrong
    (e.g., observers append the filter, and do not always update it).  Within
    each pointing group, the frames with the same filter form one mosaic.

    The labels are:

        - ``wcs_group``: the most common object name in the pointing group,
          without a trailing filter name (see :func:`_strip_filter`), made
          safe for file names (see :func:`~pipnick.metadata.names.safe_name`), e.g. ``NGC_7331``,
          with ``_1``, ``_2``, etc. appended if several pointing groups have
          the same name;
        - ``mosaic``: ``<wcs_group>_<filter>``, e.g. ``NGC_7331_V`` (see
          :func:`~pipnick.metadata.names.filter_name`).

    A filter with only one image (or combination group) in a pointing group
    is not mosaicked, so it is not labeled.  Because the labels depend only on
    the frames, the same pointing and filter get the same labels in the tables
    of different setups when they are assigned together (i.e., when the
    tables are built at the same time).

    Parameters
    ----------
    tbl : :class:`astropy.table.Table`
        The metadata table, with the ``comb_id`` column; see
        :func:`assign_comb_ids`.  The mosaic columns are added if needed (see
        :func:`add_mosaic_columns`).  Modified in place.
    commented : :class:`numpy.ndarray`, optional
        Boolean array flagging the commented-out rows, which are not assigned.
    radius : float, optional
        The linking length of the pointing groups in arcmin.  If None, half of
        the smallest field of view of the frames is used.
    """
    add_mosaic_columns(tbl)
    _commented = np.zeros(len(tbl), dtype=bool) if commented is None else commented
    unlabeled = np.array([len(str(m)) == 0 for m in tbl['mosaic']], dtype=bool)
    select = (np.asarray(tbl['frametype']) == 'science') & np.logical_not(_commented) & unlabeled
    if not np.any(select):
        return
    if radius is None:
        scale = cameras.get_camera(tbl.meta['camera']).plate_scale
        fov = np.minimum(
            np.asarray(tbl['naxis1']) * np.asarray(tbl['xbin']),
            np.asarray(tbl['naxis2']) * np.asarray(tbl['ybin'])
        )[select] * scale / 60
        _radius = 0.5 * float(np.min(fov)) * units.arcmin
    else:
        _radius = radius * units.arcmin

    # One entry per combination group, from its earliest frame
    groups = {}
    mjd = np.ma.filled(np.ma.asarray(tbl['mjd'], dtype=float), np.inf)
    for i in np.where(select)[0][np.argsort(mjd[select], kind='stable')]:
        groups.setdefault(int(tbl['comb_id'][i]), []).append(i)
    first = [rows[0] for rows in groups.values()]
    good = [
        i for i in first
        if tbl['ra'][i] is not np.ma.masked and tbl['dec'][i] is not np.ma.masked
        and len(str(tbl['ra'][i])) > 0 and len(str(tbl['dec'][i])) > 0
    ]
    if len(good) < len(first):
        log.warning(
            f'{len(first) - len(good)} science frames have no coordinates and are not assigned '
            'to mosaics.'
        )
    if len(good) == 0:
        return
    coords = SkyCoord(
        [str(tbl['ra'][i]) for i in good], [str(tbl['dec'][i]) for i in good],
        unit=(units.hourangle, units.deg)
    )
    pointing = _friends_of_friends(coords, _radius)

    names = {}
    for p in np.unique(pointing):
        objects = [
            _strip_filter(tbl['object'][i], tbl['filter'][i])
            for i, q in zip(good, pointing) if q == p and tbl['object'][i] is not np.ma.masked
        ]
        name = safe_name(max(objects, key=objects.count)) if len(objects) > 0 else 'field'
        names[p] = 'field' if len(name) == 0 else name
    counts = {n: list(names.values()).count(n) for n in set(names.values())}
    seen = {}
    labels = {}
    for p in np.unique(pointing):
        name = names[p]
        if counts[name] > 1:
            seen[name] = seen.get(name, 0) + 1
            name = f'{name}_{seen[name]}'
        labels[p] = name

    mosaic = np.array(tbl['mosaic'], dtype=object)
    wcs_group = np.array(tbl['wcs_group'], dtype=object)
    members = {}
    for i, p in zip(good, pointing):
        filt = filter_name('' if tbl['filter'][i] is np.ma.masked else tbl['filter'][i])
        members.setdefault((p, filt), []).append(int(tbl['comb_id'][i]))
    nassigned = 0
    for (p, filt), comb_ids in members.items():
        if len(comb_ids) < 2:
            continue
        rows = np.concatenate([groups[c] for c in comb_ids])
        mosaic[rows] = f'{labels[p]}_{filt}'
        wcs_group[rows] = labels[p]
        nassigned += 1
    tbl['mosaic'] = np.array(mosaic, dtype=str)
    tbl['wcs_group'] = np.array(wcs_group, dtype=str)
    log.info(
        f'Assigned {nassigned} mosaics in {len(set(labels[p] for (p, _) in members))} pointing '
        f'groups (linking length {_radius.to_value("arcmin"):.2f} arcmin).'
    )

