"""
The provenance of a reduced science image and of its sources, collected from
the products of every step.

:func:`image_record` collects what each step did with an image: its
reduction, astrometric solution, photometry, photometric calibration, and
mosaic, with the QA figures of each.  :func:`trace_source` follows one source
of an image through the source table, the removed sources, the calibrated
source table, and the standard-star measurements.  Both return a
:class:`Record`, which is rendered as text by :func:`format_record` (see
``pipnick_inspect``) and as HTML by :mod:`pipnick.report.html` (see
``pipnick_report``).  Nothing is changed.
"""

from dataclasses import dataclass, field

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.table import Table

from pipnick.datamodel.catalog import SourceCatalog
from pipnick.pkg.exceptions import PipnickDataModelError
from pipnick.pkg.paths import (
    calibrated_table, calibration_file, qa_figure, residuals_file, step_dir
)
from pipnick.report.status import (
    astrometry_status, attention_flags, calibration_status, mosaic_file, mosaic_status,
    photometry_status
)


@dataclass
class Record:
    """
    A record of the provenance of an image or a source.

    Parameters
    ----------
    title : str
        The title, e.g., the name of the image.
    sections : dict
        The sections of the record, keyed by their titles, in order.  Each is a
        list of ``(label, value)`` pairs of strings.
    figures : list
        The QA figures, as ``(caption, path)`` pairs; the figures need not
        exist.
    """
    title: str
    sections: dict = field(default_factory=dict)
    figures: list = field(default_factory=list)


def _fmt(value, fmt=None, unit=None):
    """
    Format a value for a record: ``'--'`` if it is None or masked.
    """
    if value is None or value is np.ma.masked:
        return '--'
    text = str(value) if fmt is None else format(value, fmt)
    return text if unit is None else f'{text} {unit}'


def _pm(value, error, fmt):
    """
    Format a value and its error, e.g., ``'12.345 +/- 0.012'``.
    """
    if value is None or value is np.ma.masked:
        return '--'
    if error is None or error is np.ma.masked:
        return format(value, fmt)
    return f'{format(value, fmt)} +/- {format(error, fmt)}'


def image_record(image, by_name=None, mosaic_dir=None):
    """
    Collect the provenance of a reduced science image.

    Parameters
    ----------
    image : :class:`~pipnick.datamodel.image.ReducedImage`
        The image.
    by_name : dict, optional
        The science images of the data set, keyed by file name, used to check
        a solution bootstrapped from another image; see
        :func:`~pipnick.report.status.astrometry_status`.  If None, it is not
        checked.
    mosaic_dir : str, :class:`pathlib.Path`, optional
        The directory with the mosaics; see
        :func:`~pipnick.report.status.mosaic_file`.

    Returns
    -------
    :class:`Record`
        The record, with the sections ``Reduction``, ``Astrometry``,
        ``Photometry``, ``Calibration``, and ``Mosaic``, and the QA figures of
        each step.
    """
    record = Record(image.name)
    dataset = image.dataset
    rdxdir = dataset.rdxdir
    values = image.values
    status, reason = image.reduction_status()
    reduction = [
        ('status', status if status == 'ok' else f'{status} ({reason})'),
        ('object', _fmt(values['object'])), ('filter', _fmt(values['filter'])),
        ('exposure time', _fmt(values['exptime'], '.1f', 's')),
        ('airmass', _fmt(values['airmass'], '.3f')), ('night', _fmt(image.night)),
        ('setup', dataset.tag), ('combination group', str(image.product.comb_id)),
        ('raw frames', ', '.join(f.name for f in dataset.raw_files(image.product)))
    ]
    record.sections['Reduction'] = reduction
    if status == 'missing':
        return record
    hdr = image.header
    reduction += [
        ('bias', _fmt(hdr.get('PNBIAS'))), ('flat', _fmt(hdr.get('PNFLAT'))),
        ('sky level', _fmt(hdr.get('PNSKY'), '.1f', 'e-')),
        ('cosmic-ray fraction', _fmt(hdr.get('PNCRFRAC'), '.4f')),
        ('pipnick version', _fmt(hdr.get('PNVERS')))
    ]
    record.figures.append(
        ('Reduced image', qa_figure(step_dir(rdxdir, 'reduce'), image.path.stem))
    )

    # Astrometry
    wcs = 'ok' if 'PNASTMTH' in hdr else 'missing'
    if by_name is not None:
        wcs = astrometry_status(image, by_name)
    astrometry = [('status', wcs)]
    if 'PNASTMTH' in hdr:
        astrometry += [
            ('method', hdr['PNASTMTH']), ('reference', _fmt(hdr.get('PNASTREF'))),
            ('matched sources', _fmt(hdr.get('PNASTNPT'))),
            ('rms', _fmt(hdr.get('PNASTRMS'), '.3f', 'arcsec')),
            ('pointing error', _fmt(hdr.get('PNASTPOE'), '.2f', 'arcmin'))
        ]
        if 'PNASTSUB' in hdr:
            astrometry.append(('astrometry.net submission', str(hdr['PNASTSUB'])))
    record.sections['Astrometry'] = astrometry
    record.figures.append(
        ('Astrometry', qa_figure(step_dir(rdxdir, 'astrometry'), image.path.stem))
    )

    # Photometry
    pstatus, phot = photometry_status(image)
    photometry = [('status', pstatus)]
    if phot is not None:
        psf = phot.meta.get('psf', {})
        photometry += [
            ('source table', image.phot_file.name), ('sources', str(len(phot))),
            (
                'PSF FWHM',
                f'{_fmt(psf.get("fwhm_pix"), ".2f")} pixels '
                f'({_fmt(psf.get("fwhm_arcsec"), ".2f")} arcsec)'
            ),
            (
                'PSF model',
                f'{_fmt(psf.get("fittype"))} Moffat, from {_fmt(psf.get("nstack"))} stars'
            ),
            ('aperture radius', _fmt(phot.meta.get('aper_radius_pix'), '.2f', 'pixels')),
            ('coordinates', _fmt(phot.meta.get('wcs_method')))
        ]
        removed = phot.meta.get('removed')
        if removed is not None:
            text = ', '.join(f'{n} {k}' for k, n in removed.items() if n > 0)
            total = f'{sum(removed.values())}'
            photometry.append(
                ('removed sources', total if len(text) == 0 else f'{total} ({text})')
            )
        if image.removed_file.is_file():
            photometry.append(('removed-source table', image.removed_file.name))
    record.sections['Photometry'] = photometry
    phot_dir = step_dir(rdxdir, 'photometry')
    record.figures.append(('PSF model', qa_figure(phot_dir, f'{image.path.stem}_psf')))
    record.figures.append(('Sources', qa_figure(phot_dir, f'{image.path.stem}_sources')))

    # Photometric calibration
    calibration = [('status', calibration_status(image))]
    if 'PNCALZP' in hdr:
        calibration += [
            ('file', _fmt(hdr.get('PNCALFIL'))),
            ('zero point', _pm(hdr.get('PNCALZP'), hdr.get('PNCALZE'), '.3f')),
            ('color', _fmt(hdr.get('PNCALCOL'))),
            ('color coefficient', _fmt(hdr.get('PNCALCT'), '.3f')),
            ('extinction', _fmt(hdr.get('PNCALEXT')))
        ]
    if calibrated_table(image.path).is_file():
        calibration.append(('calibrated source table', calibrated_table(image.path).name))
    if image.night is not None:
        rfile = residuals_file(rdxdir, dataset.tag, image.night)
        if rfile.is_file():
            meas = Table.read(rfile)
            meas = meas[meas['image'] == image.path.stem]
            if len(meas) > 0:
                used = int(np.sum(meas['used']))
                rejected = int(np.sum(meas['rejected']))
                calibration.append(
                    (
                        'standard measurements',
                        f'{len(meas)} ({used} used, {rejected} rejected, '
                        f'{len(meas) - used - rejected} unusable)'
                    )
                )
            record.figures.append(
                ('Calibration of the night', qa_figure(step_dir(rdxdir, 'calibrate'), rfile.stem))
            )
    record.sections['Calibration'] = calibration

    # Mosaic
    mfile = mosaic_file(image, mosaic_dir=mosaic_dir)
    mosaic = [('status', mosaic_status(image, mosaic_dir=mosaic_dir))]
    if mfile is not None:
        label, group = dataset.mosaic_labels(image.product)
        mosaic += [('mosaic', label), ('WCS group', group)]
        if mfile.is_file():
            inputs = Table.read(mfile, hdu='INPUTS')
            rows = inputs[inputs['image'] == image.name]
            if len(rows) > 0:
                row = rows[0]
                mosaic += [
                    ('reference image', 'yes' if bool(row['reference']) else 'no'),
                    ('flux scale', _fmt(row['scale'], '.4f')),
                    ('magnitude offset', _pm(row['offset'], row['offset_err'], '.4f')),
                    ('scaled', 'yes' if bool(row['scaled']) else 'no'),
                    ('PSF FWHM in', _fmt(row['fwhm_in'], '.2f', 'arcsec')),
                    ('PSF FWHM matched', _fmt(row['fwhm_out'], '.2f', 'arcsec')),
                    ('kernel', str(row['kernel'])), ('weight', _fmt(row['weight'], '.3f'))
                ]
        record.figures.append(('Mosaic', qa_figure(mfile.parent, label)))
    record.sections['Mosaic'] = mosaic
    flags = attention_flags(image, phot, dataset.camera)
    record.sections['Attention'] = [('flags', ', '.join(flags) if len(flags) > 0 else 'none')]
    return record


def trace_source(image, ra=None, dec=None, source_id=None, radius=1.):
    """
    Trace a source of an image through the products of the photometry and
    the photometric calibration.

    The source is given by its celestial coordinates, or by its ``id`` in the
    source table of the image.  A source given by its coordinates is the
    nearest source within ``radius``; if there is none, the nearest source
    removed by the photometry (see
    :func:`~pipnick.pipelines.photometry.vet_sources`) is reported instead.
    The source is then found in the calibrated source table of the image (with
    its magnitude in the other images of the same object, by its
    ``source_id``) and, for a standard star, in the standard-star measurements
    of the night, with whether it was used, rejected, or unusable, and why.

    Parameters
    ----------
    image : :class:`~pipnick.datamodel.image.ReducedImage`
        The image, which must have a source table.
    ra, dec : float, optional
        The ICRS coordinates of the source in degrees.  The image must have
        an astrometric solution.
    source_id : int, optional
        The ``id`` of the source in the source table of the image.  Used if
        ``ra`` and ``dec`` are not given.
    radius : float, optional
        The match radius in arcsec.

    Returns
    -------
    :class:`Record`
        The record, with a section for each product the source is found in.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the image has no source table, if the source is given by its
        coordinates and the image has no astrometric solution, or if neither
        the coordinates nor the ``id`` is given.
    """
    if not image.phot_file.is_file():
        raise PipnickDataModelError(f'{image.name} has no source table; run pipnick_photometry.')
    phot = SourceCatalog.read(image.phot_file)
    x = np.asarray(phot['x'], dtype=float)
    y = np.asarray(phot['y'], dtype=float)
    if ra is not None and dec is not None:
        if not image.has_solution or image.wcs is None:
            raise PipnickDataModelError(
                f'{image.name} has no astrometric solution, so a source cannot be found by its '
                'coordinates; use its id instead.'
            )
        x0, y0 = image.wcs.world_to_pixel(SkyCoord(ra, dec, unit='deg'))
        x0 = float(x0)
        y0 = float(y0)
        title = f'Source at RA {ra:.6f}, Dec {dec:.6f} in {image.name}'
    elif source_id is not None:
        rows = np.where(np.asarray(phot['id']) == int(source_id))[0]
        if len(rows) == 0:
            record = Record(f'Source {source_id} in {image.name}')
            record.sections['Not found'] = [('id', f'no source {source_id} in {phot.path.name}')]
            return record
        x0 = float(x[rows[0]])
        y0 = float(y[rows[0]])
        title = f'Source {source_id} in {image.name}'
    else:
        raise PipnickDataModelError('Give the coordinates or the id of the source.')
    record = Record(title)
    # The nominal pixel scale is accurate enough for the match radius
    scale = image.nominal_scale
    radius_pix = radius / scale
    dist = np.hypot(x - x0, y - y0)
    match = int(np.argmin(dist)) if len(dist) > 0 else None
    if match is None or dist[match] > radius_pix:
        record.sections['Source table'] = [
            ('match', f'no source within {radius:g} arcsec of pixel ({x0:.1f}, {y0:.1f})')
        ]
        _trace_removed(image, x0, y0, radius_pix, scale, record)
        return record
    row = phot[match]
    sid = int(row['id'])
    section = [
        ('table', phot.path.name), ('id', str(sid)),
        ('distance', f'{dist[match] * scale:.2f} arcsec'),
        ('pixel', f'({float(row["x"]):.2f}, {float(row["y"]):.2f})'),
        ('RA, Dec', f'{_fmt(row["ra"], ".6f")}, {_fmt(row["dec"], ".6f")}')
    ]
    # The measurements, for the columns in the table
    for label, name, fmt, unit in [
        ('PSF flux', 'flux_psf', '.1f', ' e-'), ('aperture flux', 'flux_aper', '.1f', ' e-'),
        ('instrumental magnitude', 'mag_inst', '.4f', '')
    ]:
        if name in phot.colnames:
            error = row[f'{name}_err'] if f'{name}_err' in phot.colnames else None
            section.append((label, _pm(row[name], error, fmt) + unit))
    for label, name in [('flags', 'flags'), ('detected in iteration', 'iter_detected')]:
        if name in phot.colnames:
            section.append((label, str(row[name])))
    record.sections['Source table'] = section
    _trace_calibrated(image, sid, record)
    _trace_standard(image, sid, record)
    return record


def _trace_removed(image, x0, y0, radius_pix, scale, record):
    """
    Add the nearest removed source within the radius to the record of a
    source; see :func:`trace_source`.
    """
    if not image.removed_file.is_file():
        record.sections['Removed sources'] = [
            ('table', f'{image.removed_file.name} does not exist; remeasure the image')
        ]
        return
    removed = Table.read(image.removed_file)
    dist = np.hypot(
        np.asarray(removed['x_init'], dtype=float) - x0,
        np.asarray(removed['y_init'], dtype=float) - y0
    )
    if len(dist) == 0 or np.min(dist) > radius_pix:
        record.sections['Removed sources'] = [('match', 'no removed source within the radius')]
        return
    row = removed[int(np.argmin(dist))]
    record.sections['Removed sources'] = [
        ('table', image.removed_file.name), ('reason', str(row['reason'])),
        ('distance', f'{float(np.min(dist)) * scale:.2f} arcsec'),
        ('detected at pixel', f'({float(row["x_init"]):.2f}, {float(row["y_init"]):.2f})'),
        ('fitted pixel', f'({_fmt(row["x_fit"], ".2f")}, {_fmt(row["y_fit"], ".2f")})'),
        ('fitted flux', _fmt(row['flux_fit'], '.1f', 'e-'))
    ]


def _trace_calibrated(image, sid, record):
    """
    Add a source in the calibrated source tables to its record; see
    :func:`trace_source`.
    """
    cfile = calibrated_table(image.path)
    if not cfile.is_file():
        return
    cal = Table.read(cfile)
    rows = cal[np.asarray(cal['id']) == sid]
    if len(rows) == 0:
        return
    row = rows[0]
    color = _fmt(row['color_used'], '.3f')
    if bool(row['color_assumed']):
        color += ' (assumed)'
    section = [
        ('table', cfile.name), ('band', str(cal.meta.get('band'))),
        ('magnitude', _pm(row['mag'], row['mag_err'], '.4f')),
        (f'color ({cal.meta.get("color")})', color),
        ('linked source id', str(row['source_id']))
    ]
    # The same source in the other images of the same object
    for other in sorted(cfile.parent.glob('*_cal.ecsv')):
        if other == cfile:
            continue
        tbl = Table.read(other)
        if tbl.meta.get('object') != cal.meta.get('object'):
            continue
        match = tbl[np.asarray(tbl['source_id']) == int(row['source_id'])]
        if len(match) > 0:
            section.append(
                (
                    f'in {tbl.meta.get("image")}',
                    f'{tbl.meta.get("band")} {_pm(match["mag"][0], match["mag_err"][0], ".4f")}'
                )
            )
    record.sections['Calibrated source table'] = section


def _trace_standard(image, sid, record):
    """
    Add the standard-star measurement of a source to its record; see
    :func:`trace_source`.
    """
    if image.night is None:
        return
    rfile = residuals_file(image.dataset.rdxdir, image.dataset.tag, image.night)
    if not rfile.is_file():
        return
    meas = Table.read(rfile)
    rows = meas[(meas['image'] == image.path.stem) & (np.asarray(meas['source_id']) == sid)]
    if len(rows) == 0:
        return
    row = rows[0]
    if bool(row['used']):
        use = 'used'
    elif bool(row['rejected']):
        use = 'rejected as an outlier'
    elif str(row['reason']) != '-':
        use = f'unusable: {row["reason"]}'
    else:
        use = 'not used: the band is not calibrated'
    record.sections['Standard star'] = [
        ('table', rfile.name), ('star', str(row['star'])), ('band', str(row['band'])),
        ('measurement', use),
        ('catalog magnitude', _pm(row['mag_cat'], row['mag_cat_err'], '.3f')),
        ('residual', _pm(row['resid'], row['resid_err'], '.4f')),
        ('leave-one-out residual', _fmt(row['resid_loo'], '.4f'))
    ]


def format_record(record):
    """
    Format a record as text.

    Parameters
    ----------
    record : :class:`Record`
        The record.

    Returns
    -------
    str
        The record, one ``label: value`` line per entry under the title of
        each section, followed by the QA figures and whether they exist.
    """
    lines = [record.title]
    for title, entries in record.sections.items():
        lines.append(f'  {title}:')
        width = max(len(label) for label, _ in entries)
        lines += [f'    {label:<{width}}  {value}' for label, value in entries]
    if len(record.figures) > 0:
        lines.append('  QA figures:')
        lines += [
            f'    {caption}: {path}' + ('' if path.is_file() else ' (missing)')
            for caption, path in record.figures
        ]
    return '\n'.join(lines)
