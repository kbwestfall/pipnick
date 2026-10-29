"""
The instrument setups of the frames, and the per-setup metadata files; see
:mod:`pipnick.metadata`.
"""

import operator
from pathlib import Path

import numpy as np
from astropy.table import Table

from pipnick.pkg.exceptions import PipnickDataModelError

setup_keys = ['xbin', 'ybin', 'readmode', 'window']
"""
The metadata columns that, with the camera, define an instrument setup.
"""

setup_file_keys = setup_keys + ['naxis1', 'naxis2']
"""
The metadata columns that are identical for all frames in a setup, and so are
moved from the table columns to the ``meta`` dictionary of a per-setup file;
the image size is set by the binning and readout window.
"""


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



def _window_name(index):
    """
    Return the name of a custom readout window for a 0-based index: win00,
    win01, ...
    """
    return f'win{index:02d}'



def assign_setups(tbl, camera_name, setups=None):
    """
    Assign a setup label to each frame in a metadata table.

    The ``setup`` column is (re)set, and the definition of each setup is saved
    to the ``setups`` entry of the table ``meta`` dictionary.  Each definition
    includes a ``window_name``: ``'full'`` for the full readout, or ``win00``,
    ``win01``, etc. for each distinct readout window, in order of first
    appearance.  Setups with the same window share its name.

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
        their labels, windows keep their names, and new setups and windows are
        labeled and named after them.
    """
    labels = {}
    window_names = {'full': 'full'}
    if setups is not None:
        for label, definition in setups.items():
            labels[tuple(definition[k] for k in setup_keys)] = label
            if 'window_name' in definition:
                window_names[definition['window']] = definition['window_name']
    assigned = []
    for row in tbl:
        key = tuple(row[k].item() if hasattr(row[k], 'item') else row[k] for k in setup_keys)
        if key not in labels:
            labels[key] = _setup_label(len(labels))
        window = key[setup_keys.index('window')]
        if window not in window_names:
            window_names[window] = _window_name(len(window_names) - 1)
        assigned.append(labels[key])
    tbl['setup'] = assigned
    tbl.meta['setups'] = {}
    for key, label in sorted(labels.items(), key=operator.itemgetter(1)):
        definition = {'camera': camera_name, **dict(zip(setup_keys, key))}
        definition['window_name'] = window_names[definition['window']]
        tbl.meta['setups'][label] = definition



def setup_tag(definition):
    """
    Construct the tag that identifies a setup in file names and headers.

    Parameters
    ----------
    definition : dict
        The setup definition; see :func:`assign_setups`.

    Returns
    -------
    str
        The tag, ``<xbin>x<ybin>_<readmode>_<window_name>``, e.g.,
        ``2x2_fast_full``.
    """
    return (
        f'{definition["xbin"]}x{definition["ybin"]}_{definition["readmode"]}_'
        f'{definition["window_name"]}'
    )



def setup_file_name(base, definition):
    """
    Construct the name of the metadata file for one setup.

    The name is ``<stem>_<tag><suffix>``, where ``stem`` and ``suffix`` are
    those of ``base`` and ``tag`` is from :func:`setup_tag`.

    Parameters
    ----------
    base : str, :class:`pathlib.Path`
        The base file name, e.g., ``n2.ecsv``.
    definition : dict
        The setup definition; see :func:`assign_setups`.

    Returns
    -------
    :class:`pathlib.Path`
        The file name, e.g., ``n2_2x2_fast_full.ecsv``.
    """
    _base = Path(base)
    return _base.with_name(f'{_base.stem}_{setup_tag(definition)}{_base.suffix}')



def split_by_setup(tbl, commented):
    """
    Split a metadata table into one table per setup.

    The values that are identical for all frames in a setup (see
    :data:`setup_file_keys`) and the ``setup`` column are removed from each
    table and stored in the ``setup`` entry of its ``meta`` dictionary, with the
    setup label and definition.

    Parameters
    ----------
    tbl : :class:`astropy.table.Table`
        The metadata table, with the ``setup`` column and the ``setups`` entry
        in its ``meta`` dictionary; see :func:`assign_setups`.
    commented : :class:`numpy.ndarray`
        Boolean array flagging the rows that are commented out.

    Returns
    -------
    list
        One ``(definition, table, commented)`` tuple for each setup with at
        least one frame, where ``definition`` is the ``setup`` entry of the
        table ``meta`` dictionary.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the image size is not the same for all frames in a setup.
    """
    meta = {k: v for k, v in tbl.meta.items() if k != 'setups'}
    split = []
    for label, setup in tbl.meta['setups'].items():
        indx = tbl['setup'] == label
        if not np.any(indx):
            continue
        sub = tbl[indx]
        definition = {'label': label, **setup}
        for key in ['naxis1', 'naxis2']:
            values = np.unique(sub[key])
            if len(values) != 1:
                raise PipnickDataModelError(
                    f'Frames in setup {label} have different {key} values: {list(values)}.'
                )
            definition[key] = int(values[0])
        sub.remove_columns(['setup'] + setup_file_keys)
        sub.meta = {**meta, 'setup': definition}
        split.append((definition, sub, commented[indx]))
    return split



def find_setup_files(base):
    """
    Find the existing per-setup metadata files for a base file name.

    Parameters
    ----------
    base : str, :class:`pathlib.Path`
        The base file name, e.g., ``n2.ecsv``; see :func:`setup_file_name`.

    Returns
    -------
    list
        The paths to the per-setup files, sorted by name.  Files that match the
        name pattern but are not per-setup metadata files are ignored.
    """
    _base = Path(base).absolute()
    found = []
    for f in sorted(_base.parent.glob(f'{_base.stem}_*{_base.suffix}')):
        try:
            meta = Table.read(f, format='ascii.ecsv').meta
        except (ValueError, OSError):
            continue
        if isinstance(meta.get('setup'), dict):
            found.append(f)
    return found

