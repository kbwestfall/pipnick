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
"""

import datetime
from pathlib import Path

import numpy as np
from astropy.table import Table, Column, MaskedColumn
from astropy.time import Time

from pipnick import cameras
from pipnick import log
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickPathError

__all__ = ['metadata_columns', 'setup_keys', 'build_camera_datatable', 'assign_setups',
           'build_metadata', 'frame_path', 'default_table_name', 'default_reduced_file_name']

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
}
"""
The columns of the metadata table, in order, and their data types.
"""

setup_keys = ['xbin', 'ybin', 'readmode', 'window']
"""
The metadata columns that, with the camera, define an instrument setup.
"""


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


def assign_setups(tbl, camera_name):
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
    """
    setups = {}
    labels = []
    for row in tbl:
        key = tuple(row[k].item() if hasattr(row[k], 'item') else row[k] for k in setup_keys)
        if key not in setups:
            setups[key] = _setup_label(len(setups))
        labels.append(setups[key])
    tbl['setup'] = labels
    tbl.meta['setups'] = {label: {'camera': camera_name, **dict(zip(setup_keys, key))}
                          for key, label in setups.items()}


def build_camera_datatable(root, ext='.fits', path_key='path', camera=None):
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

    Returns
    -------
    :class:`astropy.table.Table`
        Table with metadata for each frame found in ``root``, with the columns
        listed by :data:`metadata_columns`.  Rows are ordered by file name.
        The table ``meta`` dictionary includes the camera name (``camera``),
        the data path (``path_key``), and the setup definitions (``setups``).

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
    _camera = None if camera is None else cameras.get_camera(camera)

    # Get the list of files
    files = sorted(_root.glob(f'*{ext}'))
    if len(files) == 0:
        raise PipnickPathError(f'No files with extension {ext} found in {_root}.')
    log.info(f'Found {len(files)} data files in {_root}.')

    # Parse the metadata
    rows = []
    found = {}
    for f in files:
        hdr = cameras.get_header(f)
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
        rows.append({'file': f.name, **cam.metadata(hdr)})

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
    return tbl


def build_metadata(rawdir=None, filename=None, ext='.fits', overwrite=False, rdxdir=None,
                   camera=None):
    """
    Construct a table with metadata needed for the data processing.

    Execution options are:

        - If the ``filename`` is provided, it exists, and ``overwrite`` is
          false, reload the data from the file.

        - If ``rawdir`` and ``filename`` are both provided, the file provided
          exists, and overwrite is False, ``rawdir`` will be ignored and the
          execution is identical to the previous case.

        - To build the table from scratch, ``rawdir`` must be provided.  If
          ``filename`` is not provided, a default name is chosen.  If
          ``filename`` is provide and does not exist, or the file exists and
          ``overwrite`` is set to true, the provided ``filename`` is used for
          the output file.

    Parameters
    ----------
    rawdir : str, :class:`pathlib.Path`, optional
        Path to the directory containing the FITS files to be analyzed.  See
        above for the various execution modes and how this parameter interacts
        with ``filename`` and ``overwrite``.
    filename : str, optional
        The full path (if needed) and file name to either write or load with
        metadata used during the data processing.  See above for the various
        execution modes and how this parameter interacts with ``rawdir`` and
        ``overwrite``.
    ext : str, optional
        Extension of the files to include in the search string used to find
        files.
    overwrite : bool, optional
        Overwrite existing files.
    rdxdir : str, :class:`pathlib.Path`, optional
        Path to the directory for the reduced data products.  If provided, this
        is included in the table as a ``meta`` keyword.
    camera : str, optional
        The name of the camera used to obtain the data.  If None, the camera is
        identified from the file headers; see :func:`build_camera_datatable`.

    Returns
    -------
    :class:`pathlib.Path`
        Path to the file with the metadata.
    :class:`astropy.table.Table`
        Table containing file metadata.
    """
    _rdxdir = None if rdxdir is None else str(Path(rdxdir).absolute())
    _filename = None if filename is None else Path(filename).absolute()
    if _filename is None or not _filename.is_file() or _filename.is_file() and overwrite:
        # Construct the metadata table
        if rawdir is None:
            raise ValueError('Must provide a directory with the raw data if a metadata table is '
                             'not provided, it does not exist, or you wish to overwrite the file.')
        metadata = build_camera_datatable(rawdir, ext=ext, path_key='rawdir', camera=camera)
        if _filename is None:
            _filename = Path(default_table_name(metadata.meta['camera'], 'raw')).absolute()
        if _rdxdir is not None:
            metadata.meta['rdxdir'] = _rdxdir
        # Write it
        log.info(f'Saving metadata to {_filename}')
        metadata.write(_filename, format='ascii.ecsv', overwrite=True)
        # Return the table
        return _filename, metadata

    # Read and return the metadata.  This will raise an error if _filename does
    # not exist
    log.info(f'Loading metadata from {_filename}')
    metadata = Table.read(_filename, format='ascii.ecsv')
    if 'rdxdir' in metadata.meta and _rdxdir is not None and metadata.meta['rdxdir'] != _rdxdir:
        log.warning(f"Overwriting output path for reductions.  Changing from "
                    f"\n\t{metadata.meta['rdxdir']} \n to \n\t{_rdxdir}")
        metadata.meta['rdxdir'] = _rdxdir
        metadata.write(_filename, format='ascii.ecsv', overwrite=True)
    return _filename, metadata


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
