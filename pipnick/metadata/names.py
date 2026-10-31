"""
Names derived from the metadata: safe file names, filter names, and the
table names.
"""

import datetime
import re
from pathlib import Path


def safe_name(value):
    """
    Convert a string for use in a file name.

    Spaces are replaced by underscores, and any character other than letters,
    digits, underscores, periods, plus signs, and minus signs is removed.  For
    example, ``NGC 7331 I`` becomes ``NGC_7331_I``.  See :func:`filter_name`
    for filters.
    """
    return re.sub(r'[^A-Za-z0-9_.+-]', '', str(value).strip().replace(' ', '_'))



def filter_name(filt):
    """
    Convert a filter name for use in a file name.

    A prime is replaced by ``p``, so that primed filters are distinct from the
    unprimed filters with the same letter; e.g., ``i'`` becomes ``ip``, while
    ``i`` is unchanged.  The result is then converted by :func:`safe_name`, or
    is ``none`` if the filter is unknown.
    """
    name = safe_name(str(filt).replace("'", 'p'))
    return 'none' if len(name) == 0 else name



def table_file_name(name):
    """
    Construct the file name of a metadata table from the name given by the
    user.

    The name can be given with or without the ``.ecsv`` extension: ``.ecsv`` is
    appended to a name that does not end with it.  E.g., ``n3`` and
    ``n3.ecsv`` both give ``n3.ecsv``, and ``n3.night`` gives
    ``n3.night.ecsv``.  For per-setup files, the result is the base name used
    by :func:`~pipnick.metadata.setups.setup_file_name`; e.g., ``n3`` gives
    ``n3_2x2_fast_full.ecsv``.

    Parameters
    ----------
    name : str, :class:`pathlib.Path`
        The name of the table, with or without the extension.

    Returns
    -------
    :class:`pathlib.Path`
        The file name.
    """
    _name = Path(name)
    return _name if _name.suffix == '.ecsv' else _name.with_name(f'{_name.name}.ecsv')



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

