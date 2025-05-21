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
from pipnick import logger


def build_camera_datatable(root, ext='.fits', path_key='path'):
    """
    Find data files in a directory and construct a table with relevant metadata
    parsed using the relevant :class:`pipnick.cameras.cameras.NickelCamera` class.

    Parameters
    ----------
    root : str, Path
        Path to directory with the observed data files.
    ext : str, optional
        Extension of the files to use in the search string.

    Returns
    -------
    astropy.table.Table
        Table with the metadata

    Raises
    ------
    NotADirectoryError
        Raised if the provided raw directory does not exist.
    """
    # Check the inputs
    _root = Path(root).absolute()
    if not _root.is_dir():
        raise NotADirectoryError(f'{_root} does not exist!')

    # Get the list of files
    files = sorted(_root.glob(f'*{ext}'))
    nfiles = len(files)
    logger.info(f'Found {nfiles} data files.')

    # Use the first file to identify the camera
    camera_name = cameras.identify_camera(files[0])
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

#    mode=None, excl_files=None, excl_objs=None, excl_filts=None,
#    mode : str
#        Function for which organize_files() is run, & name of new table output
#        ('reduction', 'astrometry', 'photometry', 'final_calibration')
#    excl_files : list
#        List of file stems to exclude (exact match not necessary).
#    excl_objs : list
#        List of object strings to exclude (exact match not necessary).
#    excl_filts : list
#        List of filter names to exclude (exact match not necessary).

def build_metadata(rawdir=None, filename=None, ext='.fits', overwrite=False, rdxdir=None):
    """
    Construct a table with metadata needed for the data processing.

    Execution options are:

        - If the metadata has already been constructed and written to a file,
          simply pass the file to ``filename`` to reload the data.

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
    rawdir : str, Path, optional
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
    rdxdir : str, Path, optional
        Path to the directory for the reduced data products.  If provided, this
        is included in the table as a ``meta`` keyword.

    Returns
    -------
    Path
        Path to the file with the metadata.
    astropy.table.Table
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
        logger.info(f'Saving metadata to {_filename}')
        metadata.write(_filename, format='ascii.ecsv', overwrite=True)
        # Return the table
        return _filename, metadata

    # Read and return the metadata.  This will raise an error if _filename does
    # not exist
    logger.info(f'Loading metadata from {_filename}')
    metadata = Table.read(_filename, format='ascii.ecsv')
    if 'rdxdir' in metadata.meta and _rdxdir is not None and metadata.meta['rdxdir'] != _rdxdir:
        logger.warning(f"Overwriting output path for reductions.  Changing from "
                       f"\n\t{metadata.meta['rdxdir']} \n to \n\t{_rdxdir}")
        metadata.meta['rdxdir'] = _rdxdir
        metadata.write(_filename, format='ascii.ecsv', overwrite=True)
    return _filename, metadata

#    # Extract files from an astropy Table file
#    logger.info(f"Files will be extracted from Astropy table file {table_path}, not directory {datadir}")
#    # Convert astropy table to pandas DataFrame
#    file_df = file_table.to_pandas()
#    file_df.insert(1, "files", file_df.paths)
#    file_df.paths = [Path(file_path) for file_path in file_df.paths]
#    logger.info(f"{len(file_df.paths)} files extracted from table file")

#        # Create DataFrame with file metadata
#        obj_list = []
#        filt_list = []
#        for file in files:
#            hdul = fits.open(str(file))
#            obj_list.append(norm_str(hdul[0].header["OBJECT"]))
#            filt_list.append(hdul[0].header["FILTNAM"])
#            hdul.close()
#
#        file_df = pd.DataFrame({
#            "names": [file.stem for file in files],
#            "files": files,
#            "objects": obj_list,
#            "filters": filt_list,
#            "paths": files
#            })
#    
#        # Save the table file for future reference
#        logger.info(f"Saving table of {mode} file data to {table_path}")
#        logger.info(f"You can change the file name to {mode}_files.yml for ease of commenting out files in VS Code (ctrl + '/')")
#        file_table = Table.from_pandas(file_df)
#        file_table.remove_column('files')
#        file_table.write(table_path, format='ascii.fixed_width', overwrite=True)
#    
#    # Apply manual exclusions based on provided criteria
#    # (excludes file if any str in excl_list is in the the file's excl_type)
#    excl_types = ['name', 'object', 'filter', 'object']
#    excl_lists = [excl_files, excl_objs, excl_filts, [focus_label]]
#    excl_file_names_all = []
#    for excl_type, excl_list in zip(excl_types, excl_lists):
#        excl_func = create_exclusion_func(excl_list)
#        axis = file_df[excl_type+'s']
#        excl_file_names = list(file_df.names[axis.apply(lambda x: not excl_func(x))])
#        excl_file_names_all += excl_file_names
#        file_df = file_df[axis.apply(excl_func)]
#        if len(excl_file_names) > 0:
#            logger.info(f"Manually excluding files with {excl_list} in {excl_type}: {excl_file_names}")
#    
#    # Add '#' to excluded files' rows in table file to ignore them in future
#    already_excl_lines = comment_out_rows(excl_file_names_all, table_path, modify=True)
#    if len(already_excl_lines) > 0:
#        logger.info(f"Automatically excluding files already commented out in the table file: {already_excl_lines}")
#    if len(excl_file_names_all) > 0:
#        logger.info(f"Modifying table at {table_path} to ignore manually excluded files {excl_file_names_all} in future")
#    
#    # Return
#    return file_df


def frame_path(metadata, select=None, col_key='file', path_key='rawdir'):
    """
    Construct a set of file paths for rows in the provided metadata table.

    Parameters
    ----------
    metadata : astropy.table.Table
        Table with the fits file metadata
    select : numpy.ndarray, optional
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
    raw_frame : Path
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
    function
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


def unzip_directories(path_list: list, output_format: str = 'Fits_Simple', allow_exceptions: bool = False) -> list:
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
                    logger.debug(f"unzip_directories() threw KeyError or OSError on {elem_path}")
            else:
                images.append(output(elem_path))
    return images