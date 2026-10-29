"""
The summary table of a pipnick step.

Each step writes a summary of its results, with one row per item (an image, a
band, or a mosaic), to ``summary.ecsv`` in its output directory, and logs it.
Values that are unknown or do not apply are given as None, and are masked in
the table; they are never set to an invalid number.
"""

import numbers

from astropy.table import MaskedColumn, Table

from pipnick import log
from pipnick.pkg.version import version


class Summary:
    """
    The summary table of a step.

    Parameters
    ----------
    title : str
        The name of the step in the log, e.g., ``'Photometry'``.
    columns : list
        The names of the columns, in order.
    dtypes : dict, optional
        The type (e.g., ``int``, ``float``, ``str``) of each column that may
        have missing values.  Columns without missing values take the type of
        their values.
    formats : dict, optional
        The format of each column, e.g., ``'.2f'``.
    logged : list, optional
        The columns that are logged; all of them if None.
    """
    def __init__(self, title, columns, dtypes=None, formats=None, logged=None):
        self.title = title
        self.columns = list(columns)
        self.dtypes = {} if dtypes is None else dtypes
        self.formats = {} if formats is None else formats
        self.logged = self.columns if logged is None else list(logged)
        self.rows = []

    def __len__(self):
        return len(self.rows)

    def add(self, **values):
        """
        Add a row.

        Parameters
        ----------
        **values
            The values of the row, keyed by column; the columns not given are
            missing (None).

        Returns
        -------
        dict
            The row, which can be updated until the table is written.

        Raises
        ------
        KeyError
            Raised if a value is given for an unknown column.
        """
        unknown = [k for k in values if k not in self.columns]
        if len(unknown) > 0:
            raise KeyError(f'Unknown summary columns: {", ".join(unknown)}.')
        row = {name: None for name in self.columns}
        row.update(values)
        self.rows.append(row)
        return row

    def table(self):
        """
        Return the summary as a table, with the missing values masked, and the
        pipnick version in its metadata (``pipnick_version``).

        Returns
        -------
        :class:`astropy.table.Table`
            The table.
        """
        tbl = Table()
        for name in self.columns:
            values = [row[name] for row in self.rows]
            missing = [v is None for v in values]
            if any(missing) or name in self.dtypes:
                dtype = self.dtypes.get(name, _infer_type(values))
                fill = '' if dtype is str else dtype(0)
                tbl[name] = MaskedColumn(
                    [fill if v is None else dtype(v) for v in values], mask=missing, dtype=dtype
                )
            else:
                tbl[name] = values
            if name in self.formats:
                tbl[name].info.format = self.formats[name]
        tbl.meta['pipnick_version'] = version
        return tbl

    def write(self, outdir):
        """
        Write the summary to ``summary.ecsv`` in a directory, and log it.

        Parameters
        ----------
        outdir : :class:`pathlib.Path`
            The output directory of the step, created if necessary.

        Returns
        -------
        :class:`astropy.table.Table`
            The table, or None if the summary has no rows (nothing is written).
        """
        if len(self.rows) == 0:
            return None
        tbl = self.table()
        outdir.mkdir(parents=True, exist_ok=True)
        tbl.write(outdir / 'summary.ecsv', format='ascii.ecsv', overwrite=True)
        log.info(
            f'{self.title} summary (also in {outdir.name}/summary.ecsv):\n'
            + '\n'.join(tbl[self.logged].pformat(max_lines=-1, max_width=-1))
        )
        return tbl


def _infer_type(values):
    """
    Return the type of the values that are not None: ``int``, ``float``, or
    ``str``.
    """
    known = [v for v in values if v is not None]
    if len(known) > 0 and all(isinstance(v, numbers.Integral) for v in known):
        return int
    if len(known) > 0 and all(isinstance(v, numbers.Real) for v in known):
        return float
    return str
