"""
Build the metadata table of the raw frames in a directory; see
:mod:`pipnick.metadata` and :func:`build_metadata`.
"""

from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn, vstack
from astropy.time import Time

from pipnick import cameras
from pipnick import log
from pipnick.metadata.groups import add_mosaic_columns, assign_comb_ids, assign_mosaics
from pipnick.metadata.names import default_table_name, table_file_name
from pipnick.metadata.setups import (
    assign_setups, find_setup_files, setup_file_name, setup_keys, split_by_setup
)
from pipnick.metadata.table import (
    metadata_columns, read_metadata, update_metadata_meta, write_metadata
)
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickPathError
from pipnick.pkg.version import version


def _date_to_mjd(date):
    """
    Convert an ISO date string to MJD, returning None if not possible.
    """
    if date is None:
        return None
    try:
        return Time(date, format='isot', scale='utc').mjd
    except ValueError:
        return None


def _column(name, values):
    """
    Construct a typed table column, masking missing (None) values.
    """
    dtype = metadata_columns[name]
    mask = np.array([v is None for v in values], dtype=bool)
    if dtype is str:
        data = np.array(['' if v is None else v for v in values], dtype=str)
    else:
        data = np.array([0 if v is None else v for v in values], dtype=dtype)
    return MaskedColumn(data, name=name, mask=mask) if np.any(mask) else Column(data, name=name)


def build_camera_datatable(
    root, ext='.fits', path_key='path', camera=None, classify='header', skip_files=None
):
    """
    Find data files in a directory and construct a table with their metadata.

    The metadata are taken from the file headers (see
    :func:`~pipnick.cameras.cameras.NickelCamera.metadata`).  If a header does
    not provide the airmass, it is computed from the coordinates and the date
    (see :func:`~pipnick.cameras.cameras.NickelCamera.compute_airmass`), and
    ``airmass computed`` is added to the ``note`` column.

    Parameters
    ----------
    root : str, :class:`pathlib.Path`
        Path to directory with the observed data files.
    ext : str, optional
        Extension of the files to use in the search string.
    path_key : str, optional
        Metadata keyword used for the data path.
    camera : str, optional
        The name of the camera used to obtain the data (see
        :func:`~pipnick.cameras.cameras.available_cameras`).  If None, the
        camera is identified from each file header, and all files must be from
        the same camera.  If provided, files from other cameras are skipped.
    classify : str, optional
        How the frames are classified and checked.  ``'header'`` uses only the
        file headers; ``'stats'`` also uses statistics of the image data, which
        is slower but catches more problems; see
        :func:`~pipnick.cameras.cameras.NickelCamera.assess`.
    skip_files : list, optional
        Names of files in ``root`` to skip.

    Returns
    -------
    tbl : :class:`astropy.table.Table`
        Table with metadata for each frame found in ``root``, with the columns
        listed by :data:`~pipnick.metadata.table.metadata_columns`.  Rows are ordered by file name.
        The table ``meta`` dictionary includes the camera name (``camera``),
        the data path (``path_key``), the classification mode (``classify``),
        and the setup definitions (``setups``).  If no new files are found
        (all are skipped), None is returned.
    rejected : :class:`numpy.ndarray`
        Boolean array flagging the frames that failed the quality checks.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickPathError`
        Raised if the provided directory does not exist, or contains no data
        files.
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if ``camera`` is not provided and the files are from more than
        one camera, or the camera that produced a file cannot be identified.
    """
    # Check the inputs
    _root = Path(root).absolute()
    if not _root.is_dir():
        raise PipnickPathError(f'{_root} does not exist!')
    if classify not in ['header', 'stats']:
        raise ValueError(f"Classification mode must be 'header' or 'stats', not {classify}.")
    _camera = None if camera is None else cameras.get_camera(camera)

    # Get the list of files
    files = sorted(_root.glob(f'*{ext}'))
    if len(files) == 0:
        raise PipnickPathError(f'No files with extension {ext} found in {_root}.')
    if skip_files is not None:
        files = [f for f in files if f.name not in skip_files]
        if len(files) == 0:
            return None, None
    log.info(f'Found {len(files)} data files in {_root}.')

    # Parse the metadata
    rows = []
    rejected = []
    found = {}
    for f in files:
        with fits.open(f) as hdu:
            hdr = hdu[0].header
            if _camera is None:
                try:
                    cam = cameras.identify_camera(hdr)
                except PipnickDataModelError as e:
                    raise PipnickDataModelError(f'{f.name}: {e}') from e
            elif _camera.matches(hdr):
                cam = _camera
            else:
                log.warning(f'{f.name} is not from the {_camera.name} camera; skipping.')
                continue
            found.setdefault(cam.name, []).append(f.name)
            row = {'file': f.name, **cam.metadata(hdr)}
            data = hdu[0].data if classify == 'stats' else None
            row['frametype'], reject, notes = cam.assess(hdr, row['frametype'], data=data)
        if row['airmass'] is None:
            # Compute the airmass if the header does not provide it
            row['airmass'] = cam.compute_airmass(
                row['ra'], row['dec'], row['date'], exptime=row['exptime']
            )
            if row['airmass'] is not None:
                notes.append('airmass computed')
        row['note'] = '; '.join(notes)
        rows.append(row)
        rejected.append(reject)

    if len(found) == 0:
        raise PipnickPathError(f'No files from the {_camera.name} camera found in {_root}.')
    if len(found) > 1:
        summary = '; '.join(
            f'{name}: {len(names)} files (e.g., {names[0]})' for name, names in found.items()
        )
        raise PipnickDataModelError(
            f'Files in {_root} are from more than one camera '
            f'({summary}).  Select one camera, or separate the files.'
        )
    camera_name = list(found.keys())[0]
    log.info(f'Camera identified as: {camera_name}')

    # Add the MJD
    for row in rows:
        row['mjd'] = _date_to_mjd(row['date'])

    # Construct the typed table
    columns = [
        _column(name, [row.get(name) for row in rows])
        for name in metadata_columns.keys() if name != 'setup'
    ]
    tbl = Table(columns)
    assign_setups(tbl, camera_name)
    tbl = tbl[list(metadata_columns.keys())]
    tbl.meta['camera'] = camera_name
    tbl.meta[path_key] = str(_root)
    tbl.meta['classify'] = classify
    return tbl, np.array(rejected, dtype=bool)


def summarize(tbl, commented=None):
    """
    Summarize a metadata table.

    Parameters
    ----------
    tbl : :class:`astropy.table.Table`
        The metadata table.
    commented : :class:`numpy.ndarray`, optional
        Boolean array flagging the rows that are commented out.  If None, all
        rows are used.

    Returns
    -------
    str
        A summary of the setups, the number of frames of each type in each
        setup, the excluded (commented) frames, and any missing calibrations.
    """
    _commented = np.zeros(len(tbl), dtype=bool) if commented is None else commented
    use = tbl[np.logical_not(_commented)]
    frametypes = sorted(set(tbl['frametype']))
    lines = [f'Camera: {tbl.meta.get("camera")}', 'Setups:']
    for label, definition in tbl.meta.get('setups', {}).items():
        lines += [f'    {label}: ' + ', '.join(f'{k}={definition[k]}' for k in setup_keys)]

    lines += ['Frames per setup (excluding commented rows):']
    width = max(len(t) for t in frametypes + ['setup'])
    lines += ['    ' + 'setup'.ljust(6) + ' '.join(t.rjust(width) for t in frametypes)]
    for label in tbl.meta.get('setups', {}).keys():
        in_setup = use['setup'] == label
        lines += [
            '    ' + label.ljust(6) + ' '.join(
                str(np.sum(in_setup & (use['frametype'] == t))).rjust(width) for t in frametypes
            )
        ]

    nexcl = np.sum(_commented)
    lines += [f'Commented rows: {nexcl}']
    for row in tbl[_commented]:
        note = '' if row['note'] is np.ma.masked else row['note']
        lines += [f'    {row["file"]} ({row["frametype"]}): {note}']

    # Check for missing calibrations
    warnings = []
    for label in tbl.meta.get('setups', {}).keys():
        in_setup = use[use['setup'] == label]
        science = in_setup[in_setup['frametype'] == 'science']
        if len(science) == 0:
            continue
        if not np.any(in_setup['frametype'] == 'bias'):
            warnings += [f'Setup {label} has science frames but no biases.']
        flats = in_setup[np.isin(in_setup['frametype'], ['domeflat', 'skyflat'])]
        for filt in sorted(set(science['filter'])):
            if not np.any(flats['filter'] == filt):
                warnings += [f'Setup {label} has {filt} science frames but no {filt} flats.']
    if len(warnings) > 0:
        lines += ['Warnings:'] + [f'    {w}' for w in warnings]
    return '\n'.join(lines)


def _exclude_by_default(tbl, commented, select):
    """
    Comment out the pointing and focus frames, noting why.

    Parameters
    ----------
    tbl : :class:`astropy.table.Table`
        The metadata table.  The ``note`` column is modified in place.
    commented : :class:`numpy.ndarray`
        Boolean array flagging the rows that are commented out.  Modified in
        place.
    select : :class:`numpy.ndarray`
        Boolean array selecting the rows to consider, e.g., newly appended
        frames, so that earlier user edits are preserved.
    """
    exclude = select & np.isin(tbl['frametype'], ['pointing', 'focus'])
    notes = ['' if note is np.ma.masked else str(note) for note in tbl['note']]
    for i in np.where(exclude)[0]:
        reason = f'{tbl["frametype"][i]} frames are excluded by default'
        notes[i] = reason if notes[i] == '' else f'{notes[i]}; {reason}'
    tbl['note'] = notes
    commented |= exclude


def _read_setup_files(files):
    """
    Read a set of per-setup metadata files into a single table.

    Parameters
    ----------
    files : list
        The per-setup files; see :func:`~pipnick.metadata.setups.find_setup_files`.

    Returns
    -------
    tbl : :class:`astropy.table.Table`
        The combined table, with the setup columns restored, ordered by file
        name, and with the ``setups`` entry of the ``meta`` dictionary
        including all setups.
    commented : :class:`numpy.ndarray`
        Boolean array flagging the rows that are commented out.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the files are for different cameras or raw data directories.
    """
    tables = []
    commented = []
    setups = {}
    for f in files:
        tbl, com = read_metadata(f, include_commented=True)
        setups.update(tbl.meta['setups'])
        tables.append(tbl)
        commented.append(com)
    for key in ['camera', 'rawdir']:
        values = set(tbl.meta.get(key) for tbl in tables)
        if len(values) > 1:
            raise PipnickDataModelError(
                f'Per-setup files have different {key} values: {sorted(values)}.'
            )
    meta = {k: v for k, v in tables[0].meta.items() if k != 'setups'}
    combined = vstack(tables, metadata_conflicts='silent')
    commented = np.concatenate(commented)
    srt = np.argsort(combined['file'])
    combined = combined[srt]
    commented = commented[srt]
    combined.meta = {**meta, 'setups': dict(sorted(setups.items()))}
    return combined, commented


def build_metadata(
    rawdir=None, filename=None, ext='.fits', overwrite=False, append=False, rdxdir=None,
    camera=None, classify='header', split=False, assign=False, mosaic_radius=None, report=None
):
    """
    Construct a table with metadata needed for the data processing.

    Execution options are:

        - If the ``filename`` is provided, it exists, and neither ``overwrite``
          nor ``append`` is set, reload the data from the file.  ``rawdir`` is
          ignored.

        - If ``append`` is set and ``filename`` exists, files in ``rawdir``
          that are not yet in the table are added to it.  Existing rows,
          including user edits and commented-out rows, and the existing setup
          labels are preserved.

        - Otherwise, the table is built from scratch, which requires
          ``rawdir``.  If ``filename`` is not provided, a default name is
          chosen.

    If ``split`` is True, the metadata are written to one file per setup (see
    :func:`~pipnick.metadata.setups.split_by_setup` and
    :func:`~pipnick.metadata.setups.setup_file_name`), and ``filename`` is the
    base name for these files.  Existing per-setup files cannot be reloaded;
    one of ``overwrite`` or ``append`` must be set.  When appending, new frames
    are added to the file for their setup, and new setups are written to new
    files.

    Parameters
    ----------
    rawdir : str, :class:`pathlib.Path`, optional
        Path to the directory containing the FITS files to be analyzed.  See
        above for the various execution modes and how this parameter interacts
        with ``filename``, ``overwrite``, and ``append``.
    filename : str, optional
        The full path (if needed) and file name to either write or load with
        metadata used during the data processing, or the base name of the
        per-setup files if ``split`` is True.  The ``.ecsv`` extension is
        appended if the name does not end with it; see
        :func:`~pipnick.metadata.names.table_file_name`.
    ext : str, optional
        Extension of the files to include in the search string used to find
        files.
    overwrite : bool, optional
        Rebuild the table, overwriting an existing file.
    append : bool, optional
        Add new files in ``rawdir`` to an existing table.
    rdxdir : str, :class:`pathlib.Path`, optional
        Path to the directory for the reduced data products.  If provided, this
        is included in the table as a ``meta`` keyword.
    camera : str, optional
        The name of the camera used to obtain the data.  If None, the camera is
        identified from the file headers; see :func:`build_camera_datatable`.
    classify : str, optional
        How the frames are classified and checked; see
        :func:`build_camera_datatable`.
    split : bool, optional
        Write one metadata file per setup.
    assign : bool, optional
        Assign default mosaic labels to the science frames without labels;
        see :func:`~pipnick.metadata.groups.assign_mosaics`.  Only used if ``split`` is True.  When
        appending, the labels are assigned (and the files rewritten) even if
        there are no new frames.
    mosaic_radius : float, optional
        The linking length (arcmin) of the pointing groups used to assign the
        mosaics; see :func:`~pipnick.metadata.groups.assign_mosaics`.
    report : :class:`~pipnick.pkg.report.StepReport`, optional
        Records the pipnick version of an existing table that is appended
        to.  The version that writes the table is saved in its metadata
        (``pipnick_version``).

    Returns
    -------
    :class:`pathlib.Path`, list
        Path to the file with the metadata, or a list of the paths to the
        per-setup files if ``split`` is True.
    :class:`astropy.table.Table`
        Table containing file metadata for all setups, excluding commented-out
        rows.
    """
    if overwrite and append:
        raise ValueError('Cannot both overwrite and append to the metadata table.')
    _rdxdir = None if rdxdir is None else str(Path(rdxdir).absolute())
    _filename = None if filename is None else table_file_name(filename).absolute()
    if split:
        setup_files = [] if _filename is None else find_setup_files(_filename)
        exists = len(setup_files) > 0
        if exists and not overwrite and not append:
            raise ValueError(
                f'Per-setup files for {_filename} already exist; set overwrite or append.'
            )
    else:
        exists = _filename is not None and _filename.is_file()

    if exists and not overwrite and not append:
        # Read and return the metadata
        log.info(f'Loading metadata from {_filename}')
        metadata = read_metadata(_filename)
        if _rdxdir is not None and metadata.meta.get('rdxdir') != _rdxdir:
            if 'rdxdir' in metadata.meta:
                log.warning(
                    f"Overwriting output path for reductions.  Changing from "
                    f"\n\t{metadata.meta['rdxdir']} \n to \n\t{_rdxdir}"
                )
            update_metadata_meta(_filename, rdxdir=_rdxdir)
            metadata.meta['rdxdir'] = _rdxdir
        return _filename, metadata

    if rawdir is None:
        raise ValueError(
            'Must provide a directory with the raw data if a metadata table is '
            'not provided, it does not exist, or you wish to overwrite or append to '
            'the file.'
        )

    if exists and append:
        # Add new files to the existing table
        if split:
            metadata, commented = _read_setup_files(setup_files)
        else:
            metadata, commented = read_metadata(_filename, include_commented=True)
        _rawdir = str(Path(rawdir).absolute())
        if metadata.meta.get('rawdir') != _rawdir:
            raise ValueError(
                f'Cannot append files from {_rawdir} to a table of files in '
                f'{metadata.meta.get("rawdir")}.'
            )
        if metadata.meta.get('classify') != classify:
            log.warning(
                f'Appending files classified using {classify}, but existing table used '
                f'{metadata.meta.get("classify")}.'
            )
        new, rejected = build_camera_datatable(
            rawdir, ext=ext, path_key='rawdir', camera=metadata.meta['camera'], classify=classify,
            skip_files=list(metadata['file'])
        )
        if new is None and not (split and assign):
            log.info(f'No new files to add to {_filename}.')
            return setup_files if split else _filename, metadata[np.logical_not(commented)]
        if report is not None:
            report.check_version(_filename.name, metadata.meta.get('pipnick_version'))
        meta = metadata.meta.copy()
        if new is None:
            log.info(f'No new files to add to {_filename}; assigning mosaics.')
            combined = metadata
            is_new = np.zeros(len(metadata), dtype=bool)
        else:
            log.info(f'Appending {len(new)} files to {_filename}.')
            combined = vstack([metadata, new], metadata_conflicts='silent')
            commented = np.concatenate([commented, rejected])
            is_new = np.concatenate(
                [np.zeros(len(metadata), dtype=bool), np.ones(len(new), dtype=bool)]
            )
        srt = np.argsort(combined['file'])
        combined = combined[srt]
        commented = commented[srt]
        is_new = is_new[srt]
        combined.meta = meta
        assign_setups(combined, meta['camera'], setups=meta.get('setups'))
        metadata = combined
    else:
        # Construct the metadata table
        metadata, commented = build_camera_datatable(
            rawdir, ext=ext, path_key='rawdir', camera=camera, classify=classify
        )
        is_new = np.ones(len(metadata), dtype=bool)
        if _filename is None:
            _filename = Path(default_table_name(metadata.meta['camera'], 'raw')).absolute()

    if split:
        # Exclude pointing and focus frames by default, and assign the
        # combination groups of the new frames
        _exclude_by_default(metadata, commented, is_new)
        assign_comb_ids(metadata, assign=is_new)
        # Add the (empty) mosaic columns, and assign the default mosaics if
        # requested
        add_mosaic_columns(metadata)
        if assign:
            assign_mosaics(metadata, commented=commented, radius=mosaic_radius)

    if _rdxdir is not None:
        metadata.meta['rdxdir'] = _rdxdir
    metadata.meta['pipnick_version'] = version
    summary = summarize(metadata, commented)
    if not split:
        log.info(f'Saving metadata to {_filename}')
        write_metadata(metadata, _filename, commented=commented)
        log.info(f'Metadata summary:\n{summary}')
        return _filename, metadata[np.logical_not(commented)]

    # Write one file per setup
    written = []
    files_summary = ['Per-setup files:']
    for definition, tbl, com in split_by_setup(metadata, commented):
        ofile = setup_file_name(_filename, definition).absolute()
        log.info(f'Saving metadata for setup {definition["label"]} to {ofile}')
        write_metadata(tbl, ofile, commented=com)
        written.append(ofile)
        files_summary.append(
            f'    {definition["label"]}: {ofile.name} ({len(tbl)} frames, {np.sum(com)} '
            'commented)'
        )
    stale = sorted(set(setup_files) - set(written)) if exists else []
    if len(stale) > 0:
        log.warning(
            'These existing per-setup files were not rewritten and may be out of date: '
            + ', '.join(f.name for f in stale)
        )
    log.info('Metadata summary:\n' + '\n'.join([summary] + files_summary))
    return written, metadata[np.logical_not(commented)]
