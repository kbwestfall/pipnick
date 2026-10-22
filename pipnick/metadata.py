"""
Construct and manage the table of metadata for a set of observations.

The metadata table is the first product of pipnick (see the ``pipnick_metadata``
script) and defines how the data are processed.  It has one row per frame,
with the frame type, the instrument setup, and the header quantities needed
for the reduction, normalized across the Nickel cameras (see
:mod:`pipnick.cameras`).  The table is written as an ECSV file.

Setups
------

Frames obtained with the same camera, binning, readout mode, and readout
window share a *setup*, labeled ``A``, ``B``, etc. in order of first
appearance.  A typical night uses a single setup.  Calibrations (biases,
darks, and flats) can only be applied to frames with the same setup.  The
definition of each setup is stored in the ``setups`` entry of the table
``meta`` dictionary.

Excluding frames
----------------

Frames are excluded from the reduction by *commenting out* their rows in the
ECSV file, i.e., by adding a ``#`` to the beginning of the line.  Frames that
fail the quality checks (see
:func:`~pipnick.cameras.cameras.NickelCamera.assess`) are commented out
automatically, with the reason given in the ``note`` column.  Users can
comment out (or restore) any row, and can edit the table, e.g., to correct a
frame type.  These edits are preserved when new frames are appended to the
table; see :func:`build_metadata`.
"""

import csv
import datetime
import io
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn, vstack
from astropy.time import Time

from pipnick import cameras
from pipnick import log
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickPathError

__all__ = ['metadata_columns', 'setup_keys', 'classify_modes', 'build_camera_datatable',
           'assign_setups', 'write_metadata', 'read_metadata', 'update_metadata_meta',
           'summarize', 'build_metadata', 'frame_path', 'default_table_name',
           'default_reduced_file_name']

metadata_columns = {
    'file': str,
    'frametype': str,
    'setup': str,
    'object': str,
    'filter': str,
    'exptime': float,
    'airmass': float,
    'ra': str,
    'dec': str,
    'date': str,
    'mjd': float,
    'xbin': int,
    'ybin': int,
    'readmode': str,
    'naxis1': int,
    'naxis2': int,
    'window': str,
    'note': str,
}
"""
The columns of the metadata table, in order, and their data types.
"""

setup_keys = ['xbin', 'ybin', 'readmode', 'window']
"""
The metadata columns that, with the camera, define an instrument setup.
"""

classify_modes = ['header', 'stats']
"""
How frames are classified and checked.  ``'header'`` uses only the file
headers; ``'stats'`` also uses statistics of the image data, which is slower
but catches more problems; see
:func:`~pipnick.cameras.cameras.NickelCamera.assess`.
"""


def _split_ecsv_line(line):
    """
    Split a line of ECSV data into its values (space-delimited, with double
    quotes around values that contain spaces).
    """
    return next(csv.reader([line], delimiter=' ', quotechar='"'))


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


def _setup_label(index):
    """
    Return the setup label for a 0-based index: A, B, ..., Z, AA, AB, ...
    """
    label = ''
    index += 1
    while index > 0:
        index, rem = divmod(index - 1, 26)
        label = chr(ord('A') + rem) + label
    return label


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


def assign_setups(tbl, camera_name, setups=None):
    """
    Assign a setup label to each frame in a metadata table.

    The ``setup`` column is (re)set, and the definition of each setup is saved
    to the ``setups`` entry of the table ``meta`` dictionary.

    Parameters
    ----------
    tbl : :class:`astropy.table.Table`
        The metadata table.  Must contain the :data:`setup_keys` columns.
        Modified in place.
    camera_name : str
        The name of the camera used to obtain the data.
    setups : dict, optional
        Existing setup definitions, as stored in the ``setups`` entry of a
        metadata table ``meta`` dictionary.  Frames with these setups keep
        their labels, and new setups are labeled after them.
    """
    labels = {}
    if setups is not None:
        for label, definition in setups.items():
            labels[tuple(definition[k] for k in setup_keys)] = label
    assigned = []
    for row in tbl:
        key = tuple(row[k].item() if hasattr(row[k], 'item') else row[k] for k in setup_keys)
        if key not in labels:
            labels[key] = _setup_label(len(labels))
        assigned.append(labels[key])
    tbl['setup'] = assigned
    tbl.meta['setups'] = {label: {'camera': camera_name, **dict(zip(setup_keys, key))}
                          for key, label in sorted(labels.items(), key=lambda x: x[1])}


def build_camera_datatable(root, ext='.fits', path_key='path', camera=None, classify='header',
                           skip_files=None):
    """
    Find data files in a directory and construct a table with their metadata.

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
        How the frames are classified and checked; see :data:`classify_modes`.
    skip_files : list, optional
        Names of files in ``root`` to skip.

    Returns
    -------
    tbl : :class:`astropy.table.Table`
        Table with metadata for each frame found in ``root``, with the columns
        listed by :data:`metadata_columns`.  Rows are ordered by file name.
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
    if classify not in classify_modes:
        raise ValueError(f'Classification mode must be one of {classify_modes}, not {classify}.')
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
        row['note'] = '; '.join(notes)
        rows.append(row)
        rejected.append(reject)

    if len(found) == 0:
        raise PipnickPathError(f'No files from the {_camera.name} camera found in {_root}.')
    if len(found) > 1:
        summary = '; '.join(f'{name}: {len(names)} files (e.g., {names[0]})'
                            for name, names in found.items())
        raise PipnickDataModelError(f'Files in {_root} are from more than one camera '
                                    f'({summary}).  Select one camera, or separate the files.')
    camera_name = list(found.keys())[0]
    log.info(f'Camera identified as: {camera_name}')

    # Add the MJD
    for row in rows:
        row['mjd'] = _date_to_mjd(row['date'])

    # Construct the typed table
    tbl = Table([_column(name, [row.get(name) for row in rows])
                 for name in metadata_columns.keys() if name != 'setup'])
    assign_setups(tbl, camera_name)
    tbl = tbl[list(metadata_columns.keys())]
    tbl.meta['camera'] = camera_name
    tbl.meta[path_key] = str(_root)
    tbl.meta['classify'] = classify
    return tbl, np.array(rejected, dtype=bool)


def write_metadata(tbl, filename, commented=None, overwrite=True):
    """
    Write a metadata table to an ECSV file, commenting out selected rows.

    Parameters
    ----------
    tbl : :class:`astropy.table.Table`
        The metadata table.
    filename : str, :class:`pathlib.Path`
        The output file.
    commented : :class:`numpy.ndarray`, optional
        Boolean array selecting the rows to comment out.  If None, no rows are
        commented out.
    overwrite : bool, optional
        Overwrite an existing file.
    """
    _filename = Path(filename).absolute()
    if _filename.exists() and not overwrite:
        raise FileExistsError(f'{_filename} exists and overwrite is False.')
    buf = io.StringIO()
    tbl.write(buf, format='ascii.ecsv')
    lines = buf.getvalue().splitlines()
    if commented is not None and np.any(commented):
        # The data rows follow the line with the column names, which is the
        # first line that is not part of the header
        start = next(i for i, line in enumerate(lines) if not line.startswith('#')) + 1
        for i in np.where(commented)[0]:
            lines[start + i] = '# ' + lines[start + i]
    _filename.write_text('\n'.join(lines) + '\n')


def read_metadata(filename, include_commented=False):
    """
    Read a metadata table from an ECSV file.

    Parameters
    ----------
    filename : str, :class:`pathlib.Path`
        The ECSV file.
    include_commented : bool, optional
        Include the commented-out rows in the returned table.  Comment lines
        that are not table rows are ignored.

    Returns
    -------
    tbl : :class:`astropy.table.Table`
        The metadata table.  If ``include_commented`` is False, the
        commented-out rows are excluded.
    commented : :class:`numpy.ndarray`
        Boolean array flagging the rows that were commented out.  Only returned
        if ``include_commented`` is True.
    """
    _filename = Path(filename).absolute()
    if not include_commented:
        return Table.read(_filename, format='ascii.ecsv')

    lines = _filename.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if not line.startswith('#'))
    ncol = len(_split_ecsv_line(lines[start]))
    commented = []
    for i in range(start + 1, len(lines)):
        if not lines[i].startswith('#'):
            commented.append(False)
            continue
        candidate = lines[i].lstrip('#').lstrip()
        if len(_split_ecsv_line(candidate)) == ncol:
            lines[i] = candidate
            commented.append(True)
    tbl = Table.read('\n'.join(lines), format='ascii.ecsv')
    return tbl, np.array(commented, dtype=bool)


def update_metadata_meta(filename, **kwargs):
    """
    Update the ``meta`` dictionary of a metadata table file.

    Commented-out rows are preserved.

    Parameters
    ----------
    filename : str, :class:`pathlib.Path`
        The ECSV file, which is overwritten.
    **kwargs
        The ``meta`` entries to set.
    """
    tbl, commented = read_metadata(filename, include_commented=True)
    tbl.meta.update(kwargs)
    write_metadata(tbl, filename, commented=commented)


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
        lines += ['    ' + label.ljust(6)
                  + ' '.join(str(np.sum(in_setup & (use['frametype'] == t))).rjust(width)
                             for t in frametypes)]

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


def build_metadata(rawdir=None, filename=None, ext='.fits', overwrite=False, append=False,
                   rdxdir=None, camera=None, classify='header'):
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

    Parameters
    ----------
    rawdir : str, :class:`pathlib.Path`, optional
        Path to the directory containing the FITS files to be analyzed.  See
        above for the various execution modes and how this parameter interacts
        with ``filename``, ``overwrite``, and ``append``.
    filename : str, optional
        The full path (if needed) and file name to either write or load with
        metadata used during the data processing.
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
        How the frames are classified and checked; see :data:`classify_modes`.

    Returns
    -------
    :class:`pathlib.Path`
        Path to the file with the metadata.
    :class:`astropy.table.Table`
        Table containing file metadata, excluding commented-out rows.
    """
    if overwrite and append:
        raise ValueError('Cannot both overwrite and append to the metadata table.')
    _rdxdir = None if rdxdir is None else str(Path(rdxdir).absolute())
    _filename = None if filename is None else Path(filename).absolute()
    exists = _filename is not None and _filename.is_file()

    if exists and not overwrite and not append:
        # Read and return the metadata
        log.info(f'Loading metadata from {_filename}')
        metadata = read_metadata(_filename)
        if _rdxdir is not None and metadata.meta.get('rdxdir') != _rdxdir:
            if 'rdxdir' in metadata.meta:
                log.warning(f"Overwriting output path for reductions.  Changing from "
                            f"\n\t{metadata.meta['rdxdir']} \n to \n\t{_rdxdir}")
            update_metadata_meta(_filename, rdxdir=_rdxdir)
            metadata.meta['rdxdir'] = _rdxdir
        return _filename, metadata

    if rawdir is None:
        raise ValueError('Must provide a directory with the raw data if a metadata table is '
                         'not provided, it does not exist, or you wish to overwrite or append to '
                         'the file.')

    if exists and append:
        # Add new files to the existing table
        metadata, commented = read_metadata(_filename, include_commented=True)
        _rawdir = str(Path(rawdir).absolute())
        if metadata.meta.get('rawdir') != _rawdir:
            raise ValueError(f'Cannot append files from {_rawdir} to a table of files in '
                             f'{metadata.meta.get("rawdir")}.')
        if metadata.meta.get('classify') != classify:
            log.warning(f'Appending files classified using {classify}, but existing table used '
                        f'{metadata.meta.get("classify")}.')
        new, rejected = build_camera_datatable(rawdir, ext=ext, path_key='rawdir',
                                               camera=metadata.meta['camera'], classify=classify,
                                               skip_files=list(metadata['file']))
        if new is None:
            log.info(f'No new files to add to {_filename}.')
            return _filename, metadata[np.logical_not(commented)]
        log.info(f'Appending {len(new)} files to {_filename}.')
        meta = metadata.meta.copy()
        combined = vstack([metadata, new], metadata_conflicts='silent')
        commented = np.concatenate([commented, rejected])
        srt = np.argsort(combined['file'])
        combined = combined[srt]
        commented = commented[srt]
        combined.meta = meta
        assign_setups(combined, meta['camera'], setups=meta.get('setups'))
        metadata = combined
    else:
        # Construct the metadata table
        metadata, commented = build_camera_datatable(rawdir, ext=ext, path_key='rawdir',
                                                     camera=camera, classify=classify)
        if _filename is None:
            _filename = Path(default_table_name(metadata.meta['camera'], 'raw')).absolute()

    if _rdxdir is not None:
        metadata.meta['rdxdir'] = _rdxdir
    log.info(f'Saving metadata to {_filename}')
    write_metadata(metadata, _filename, commented=commented)
    log.info(f'Metadata summary:\n{summarize(metadata, commented)}')
    return _filename, metadata[np.logical_not(commented)]


def frame_path(metadata, select=None, col_key='file', path_key='rawdir'):
    """
    Construct a set of file paths for rows in the provided metadata table.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        Table with the fits file metadata
    select : :class:`numpy.ndarray`, optional
        A boolean or integer vector selecting the rows of the table to use.  If
        None, all rows are used.
    col_key : str, optional
        The keyword of the table column with the file names.
    path_key : str, optional
        The keyword in the table ``meta`` dictionary with the root path of each
        file.  Cannot be None, and must be in the ``meta`` dictionary of the
        provided table.

    Returns
    -------
    list
        A list of Path objects with the full path to each file.
    """
    if path_key not in metadata.meta:
        raise ValueError(f'Provided table does not include {path_key} in its metadata dictionary.')
    files = metadata[col_key] if select is None else metadata[col_key][select]
    return [Path(metadata.meta[path_key]).absolute() / f for f in files]


def default_table_name(camera, name):
    """
    Construct the default name for an output datatable.

    The name of the file is ``{camera}_{time}_{name}.tbl``.

    Parameters
    ----------
    camera : str
        The name of the camera used to take the data.
    name : str
        A generic signifier for the type of data table.  E.g., ``'raw'`` for the
        table with the metadata for the raw observations.

    Returns
    -------
    str
        The file name.
    """
    dtime = datetime.datetime.now(datetime.UTC).strftime('%Y-%m-%dT%H:%M:%S')
    return f'{camera}_{dtime}_{name}.tbl'


def default_reduced_file_name(raw_frame):
    """
    Return the default name of a reduced image based on the name of raw file.

    Parameters
    ----------
    raw_frame : :class:`pathlib.Path`
        Path to the raw file.

    Returns
    -------
    str
        Name for the reduced data file.
    """
    return f'{raw_frame.stem}_rdx.fits'
