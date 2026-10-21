from pathlib import Path
import datetime

from IPython import embed


import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

from pipnick.utils.fits_class import Fits_Simple
from pipnick.utils.nickel_data import focus_label
from pipnick import cameras
from pipnick import log


def build_camera_datatable(root, ext='.fits', path_key='path'):
    """
    Find data files in a directory and construct a table with relevant metadata
    parsed using the relevant :class:`~pipnick.cameras.cameras.NickelCamera`
    class.

    Parameters
    ----------
    root : str, :class:`pathlib.Path`
        Path to directory with the observed data files.
    ext : str, optional
        Extension of the files to use in the search string.
    path_key : str, optional
        Metadata keyword used for the data path.

    Returns
    -------
    :class:`astropy.table.Table`
        Table with metadata for each frame found in ``root``.  The columns are
        defined by the ``parse_metadata`` function of the class relevant to
        each camera.  The camera is identified by the header of the first
        (alphabetically) raw file found in ``root``.

    Raises
    ------
    :exc:`NotADirectoryError`
        Raised if the provided raw directory does not exist.
    """
    # Check the inputs
    _root = Path(root).absolute()
    if not _root.is_dir():
        raise NotADirectoryError(f'{_root} does not exist!')

    # Get the list of files
    files = sorted(_root.glob(f'*{ext}'))
    nfiles = len(files)
    log.info(f'Found {nfiles} data files.')

    # Use the first file to identify the camera
    log.info(f'Using {files[0].name} to identify the camera.')
    camera_name = cameras.identify_camera(files[0])
    log.info(f'Camera identified as: {camera_name}')
    camera = eval(f'cameras.{camera_name}')

    # Parse the metadata into a table and write the file
    metadata = [None]*nfiles
    cols = camera.metadata_cols()
    for i,f in enumerate(files):
        metadata[i] = camera.parse_metadata(f)
    tbl = Table(data=np.asarray(metadata), names=cols)
    tbl.meta['camera'] = camera_name
    tbl.meta[path_key] = str(_root)
    return tbl


def build_metadata(rawdir=None, filename=None, ext='.fits', overwrite=False, rdxdir=None):
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
        metadata = build_camera_datatable(rawdir, ext=ext, path_key='rawdir')
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


def comment_out_rows(excluded_file_names, table_file, modify=True):
    """
    Comment out specified rows in a table file based on exclusion criteria.

    Parameters
    ----------
    excluded_file_names : list
        List of file names to comment out.
    table_file : str
        Path to the table file to modify.
    modify : bool, optional
        Whether to modify the file or not.

    Returns
    -------
    list
        List of file names that were already commented out.
    """
    with open(table_file, 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    already_excl_lines = []
    for line in lines:
        if line.strip().startswith('#'):
            # Retrieves just the file stem (i.e. 'd1001')
            already_excl_lines.append(line.split('|')[1].split(' ')[1])
            new_lines.append(line)
        elif modify and any(file_name in line for file_name in excluded_file_names):
            new_lines.append('#' + line)    # "comments out" the row
        elif modify:
            new_lines.append(line)
    
    if modify:
        with open(table_file, 'w') as f:
            f.writelines(new_lines)
    
    return already_excl_lines


def norm_str(s):
    """
    Normalize a string for comparison purposes--all caps, no spaces.
    'Sky flat' -> 'SKYFLAT'

    Parameters
    ----------
    s : str or list
        String or list of strings to normalize.

    Returns
    -------
    str or list
        Normalized string or list of normalized strings.
    """
    if isinstance(s, list):
        return [norm_str(elem) for elem in s]
    return s.upper().replace(' ', '')

def create_exclusion_func(exclude_list):
    """
    Create a function to determine if any string in the provided
    exclude_list is in another string, for the purpose of excluding files.
    
    Parameters
    ----------
    exclude_list : list
        List of strings for exclusion.

    Returns
    -------
    callable
        Function that takes a target (string) and returns True if any
        string in exclude_list is in target.
    """
    if exclude_list is None:
        return lambda _: True
    exclude_list = [norm_str(obj_str) for obj_str in exclude_list]
    def excl_func(target):
        target_str = norm_str(target)
        is_excluded = any(excluded_str in target_str for excluded_str in exclude_list)
        return not is_excluded
    return excl_func


def unzip_directories(path_list: list, output_format: str = 'Fits_Simple',
                      allow_exceptions: bool = False) -> list:
    """
    Extract all files in a list of directories or files (or their subdirectories) and
    return a list of image objects.

    Parameters
    ----------
    path_list : list
        A list of paths to directories or files that need to be unzipped.
    output_format : str, optional
        The format of the output objects. Can be 'Fits_Simple' (default) or 'Path'.
    allow_exceptions : bool, optional
        Whether to allow and handle exceptions for file-related errors such
        as `KeyError` or `OSError`. Defaults to False.

    Returns
    -------
    list
        A list of image objects in the specified output format. The objects can be
        instances of `Fits_Simple` or `Path`, depending on input `output_format`.
    """
    if output_format == 'Path':
        output = Path
    elif output_format == 'Fits_Simple':
        output = Fits_Simple
    
    images = []
    for elem_path in path_list:
        elem_path = Path(elem_path)
        if elem_path.is_dir():
            for sub_elem_path in elem_path.iterdir():
                images += unzip_directories([sub_elem_path], output_format, allow_exceptions)
        elif elem_path.is_file():
            if allow_exceptions:
                try:
                    images.append(output(elem_path))
                except (KeyError, OSError):
                    log.debug(f"unzip_directories() threw KeyError or OSError on {elem_path}")
            else:
                images.append(output(elem_path))
    return images
