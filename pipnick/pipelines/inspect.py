"""
Inspection of a reduced science image and its sources.

:func:`inspect_image` logs the provenance of one image (see
:func:`~pipnick.report.records.image_record`), optionally traces one of its
sources through the photometry and calibration products (see
:func:`~pipnick.report.records.trace_source`), and optionally shows the
image with its sources.  It only reads the products.
"""

from pathlib import Path

import numpy as np

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.table import Table

from pipnick import log
from pipnick.datamodel.catalog import SourceCatalog
from pipnick.datamodel.dataset import DataSet
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickParameterError
from pipnick.qa.inspect import plot_inspection
from pipnick.report.records import format_record, image_record, trace_source


def parse_coordinates(text):
    """
    Parse the coordinates of a source given on the command line.

    Parameters
    ----------
    text : str
        ``<RA>,<Dec>``, in degrees (e.g., ``282.77,-6.27``) or sexagesimal,
        with the RA in hours (e.g., ``18:51:05.0,-06:16:12``).

    Returns
    -------
    ra, dec : float
        The coordinates in degrees.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if the coordinates cannot be parsed.
    """
    parts = [p.strip() for p in str(text).split(',')]
    if len(parts) != 2:
        raise PipnickParameterError(f'Invalid coordinates {text}; must be RA,Dec.')
    try:
        if ':' in parts[0]:
            coo = SkyCoord(parts[0], parts[1], unit=(units.hourangle, units.deg))
        else:
            coo = SkyCoord(float(parts[0]), float(parts[1]), unit='deg')
    except ValueError as e:
        raise PipnickParameterError(f'Invalid coordinates {text}: {e}') from e
    return float(coo.ra.deg), float(coo.dec.deg)


def inspect_image(
    table, image, rdxdir=None, mosaic_dir=None, ra=None, dec=None, source_id=None, radius=1.,
    display=False, ofile=None
):
    """
    Log the provenance of a reduced science image and, optionally, of one of
    its sources, and show the image with its sources.

    Parameters
    ----------
    table : str, :class:`pathlib.Path`
        A per-setup metadata file.
    image : str
        The image: its file name, its file name without the extension, or its
        combination group (``comb_id``).
    rdxdir : str, :class:`pathlib.Path`, optional
        The reduction directory.  If None, the directory saved in the metadata
        file is used.
    mosaic_dir : str, :class:`pathlib.Path`, optional
        The directory with the mosaics.  If None, the ``mosaic`` subdirectory
        of the reduction directory is used.
    ra, dec : float, optional
        The ICRS coordinates in degrees of a source to trace; see
        :func:`~pipnick.report.records.trace_source`.
    source_id : int, optional
        The ``id`` of a source to trace, in the source table of the image.
        Used if ``ra`` and ``dec`` are not given.
    radius : float, optional
        The match radius in arcsec for a source given by its coordinates.
    display : bool, optional
        Show the image with its sources (see
        :func:`~pipnick.qa.inspect.plot_inspection`) in an interactive window.
    ofile : str, :class:`pathlib.Path`, optional
        Write the image with its sources to this file instead of showing it.

    Returns
    -------
    record : :class:`~pipnick.report.records.Record`
        The provenance of the image.
    trace : :class:`~pipnick.report.records.Record`
        The provenance of the source, or None if no source is traced.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if ``image`` does not select exactly one science image.
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the image must be shown, or a source traced, but it has not
        been reduced or has no source table.
    """
    dataset = DataSet.from_file(table, rdxdir=rdxdir, read_only=True)
    try:
        selected = dataset.science_images(selection=[image], include_missing=True)
    except ValueError as e:
        raise PipnickParameterError(str(e)) from e
    if len(selected) != 1:
        raise PipnickParameterError(
            f'{image} selects {len(selected)} images; give the file name of one image.'
        )
    target = selected[0]
    by_name = {i.name: i for i in dataset.science_images(include_missing=True)}
    record = image_record(target, by_name=by_name, mosaic_dir=mosaic_dir)
    log.info(format_record(record))

    trace = None
    if (ra is not None and dec is not None) or source_id is not None:
        trace = trace_source(target, ra=ra, dec=dec, source_id=source_id, radius=radius)
        log.info(format_record(trace))
    if not display and ofile is None:
        return record, trace

    if not target.exists:
        raise PipnickDataModelError(f'{target.name} has not been reduced; run pipnick_reduce.')
    if not target.phot_file.is_file():
        raise PipnickDataModelError(f'{target.name} has no source table; run pipnick_photometry.')
    data, mask, _, _ = target.read()
    phot = SourceCatalog.read(target.phot_file)
    removed = Table.read(target.removed_file) if target.removed_file.is_file() else None
    highlight = None
    if ra is not None and dec is not None and target.has_solution:
        x, y = target.wcs.world_to_pixel(SkyCoord(ra, dec, unit='deg'))
        highlight = (float(x), float(y))
    elif source_id is not None:
        rows = np.where(np.asarray(phot['id']) == int(source_id))[0]
        if len(rows) > 0:
            highlight = (float(phot['x'][rows[0]]), float(phot['y'][rows[0]]))
    title = f'{target.name}: {len(phot)} sources, labeled by id'
    plot_inspection(
        data, mask, phot.table, removed=removed, highlight=highlight, title=title,
        ofile=None if ofile is None else Path(ofile)
    )
    if ofile is not None:
        log.info(f'Wrote {ofile}')
    return record, trace
