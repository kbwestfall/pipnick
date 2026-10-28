"""
Read and write the metadata table, preserving its commented-out rows; see
:mod:`pipnick.metadata`.
"""

import csv
import io
from pathlib import Path

import numpy as np
from astropy.table import Column, Table

from pipnick.metadata.groups import add_mosaic_columns
from pipnick.metadata.setups import setup_file_keys

metadata_columns = {
    'file': str, 'frametype': str, 'setup': str, 'object': str, 'filter': str, 'exptime': float,
    'airmass': float, 'ra': str, 'dec': str, 'date': str, 'mjd': float, 'xbin': int,
    'ybin': int, 'readmode': str, 'naxis1': int, 'naxis2': int, 'window': str, 'note': str
}
"""
The columns of the metadata table, in order, and their data types.
"""


def _split_ecsv_line(line):
    """
    Split a line of ECSV data into its values (space-delimited, with double
    quotes around values that contain spaces).
    """
    return next(csv.reader([line], delimiter=' ', quotechar='"'))



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



def read_metadata(filename, include_commented=False, expand=True):
    """
    Read a metadata table from an ECSV file.

    Parameters
    ----------
    filename : str, :class:`pathlib.Path`
        The ECSV file, with either all setups or a single setup; see
        :func:`~pipnick.metadata.setups.split_by_setup`.
    include_commented : bool, optional
        Include the commented-out rows in the returned table.  Comment lines
        that are not table rows are ignored.
    expand : bool, optional
        For a per-setup file, restore the setup columns from its ``meta``
        dictionary, so that the table has the same layout as a table with all
        setups.  If False, the table is returned as it is in the file.

    Returns
    -------
    tbl : :class:`astropy.table.Table`
        The metadata table.  If ``include_commented`` is False, the
        commented-out rows are excluded.  A per-setup file always has the
        mosaic columns (see :data:`~pipnick.metadata.groups.mosaic_columns`),
        with empty strings for empty labels.
    commented : :class:`numpy.ndarray`
        Boolean array flagging the rows that were commented out.  Only returned
        if ``include_commented`` is True.
    """
    _filename = Path(filename).absolute()
    if not include_commented:
        return _read_layout(Table.read(_filename, format='ascii.ecsv'), expand)

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
    return _read_layout(tbl, expand), np.array(commented, dtype=bool)



def _read_layout(tbl, expand):
    """
    Finish reading a metadata table: restore the setup columns of a per-setup
    file if ``expand`` is True (see :func:`_expand_setup`), and add or fill its
    mosaic columns (see :func:`~pipnick.metadata.groups.add_mosaic_columns`),
    whose empty labels are read as masked values.
    """
    _tbl = _expand_setup(tbl) if expand else tbl
    if 'comb_id' in _tbl.colnames:
        add_mosaic_columns(_tbl)
    return _tbl



def _expand_setup(tbl):
    """
    Restore the setup columns of a table read from a per-setup file.

    The values in the ``setup`` entry of the ``meta`` dictionary are added as
    columns, in their usual positions; the entry is replaced by a ``setups``
    entry with the single setup definition, as in a table with all setups.
    Tables without a ``setup`` entry are returned unchanged.
    """
    if not isinstance(tbl.meta.get('setup'), dict):
        return tbl
    definition = tbl.meta.pop('setup')
    values = {'setup': definition['label'], **{k: definition[k] for k in setup_file_keys}}
    for name in metadata_columns.keys():
        if name not in values or name in tbl.colnames:
            continue
        # Insert after the closest preceding column that is already present
        preceding = list(metadata_columns.keys())[:list(metadata_columns.keys()).index(name)]
        present = [c for c in preceding if c in tbl.colnames]
        index = 0 if len(present) == 0 else tbl.colnames.index(present[-1]) + 1
        dtype = metadata_columns[name]
        # NOTE: Setting dtype=str would truncate the strings to one character, so
        # the string length is instead inferred from the value
        data = (
            np.full(len(tbl), values[name]) if dtype is str
            else np.full(len(tbl), values[name], dtype=dtype)
        )
        tbl.add_column(Column(data, name=name), index=index)
    setup = {k: v for k, v in definition.items() if k not in ['label', 'naxis1', 'naxis2']}
    tbl.meta['setups'] = {definition['label']: setup}
    return tbl



def update_metadata_meta(filename, **kwargs):
    """
    Update the ``meta`` dictionary of a metadata table file.

    Commented-out rows are preserved, as is the layout of a per-setup file.

    Parameters
    ----------
    filename : str, :class:`pathlib.Path`
        The ECSV file, which is overwritten.
    **kwargs
        The ``meta`` entries to set.
    """
    tbl, commented = read_metadata(filename, include_commented=True, expand=False)
    tbl.meta.update(kwargs)
    write_metadata(tbl, filename, commented=commented)



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

