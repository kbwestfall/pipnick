"""
Report the status of the processing of a data set; see
:mod:`pipnick.report.status`.
"""

from pipnick.report.status import log_status, night_status


def status_all(table, rdxdir=None, mosaic_dir=None, ofile=None):
    """
    Log the status of the processing of the science images of a data set,
    and the next command to run.

    Nothing is changed; see :func:`~pipnick.report.status.night_status`.

    Parameters
    ----------
    table : str, :class:`pathlib.Path`
        A per-setup metadata file.
    rdxdir : str, :class:`pathlib.Path`, optional
        The reduction directory.  If None, the directory saved in the metadata
        file is used.
    mosaic_dir : str, :class:`pathlib.Path`, optional
        The directory with the mosaics.  If None, the ``mosaic`` subdirectory
        of the reduction directory is used.
    ofile : str, :class:`pathlib.Path`, optional
        Also write the status table to this ECSV file.

    Returns
    -------
    tbl : :class:`astropy.table.Table`
        The status table; see :func:`~pipnick.report.status.night_status`.
    next_step : str
        The next command to run.
    """
    tbl, next_step = night_status(table, rdxdir=rdxdir, mosaic_dir=mosaic_dir)
    log_status(tbl, next_step)
    if ofile is not None:
        tbl.write(ofile, format='ascii.ecsv', overwrite=True)
    return tbl, next_step
