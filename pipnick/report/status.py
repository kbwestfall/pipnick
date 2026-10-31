"""
The status of the processing of a data set.

:func:`night_status` reports, for each science image of a per-setup metadata
table, the state of each processing step, from the products and their
provenance; nothing is changed.  Each step has one of these states:

    - ``ok``: the product exists and is up to date;
    - ``stale``: the product exists, but its inputs have changed, so the
      step must be run again;
    - ``missing``: the product does not exist;
    - ``failed``: the step could not construct the product (the photometry);
    - ``n/a``: the step does not apply, e.g., the image is not in a mosaic,
      or the night has no calibration.

The status of a step is ``n/a`` if a step it depends on is ``missing``.
"""

from pathlib import Path

import numpy as np

from astropy.io import fits
from astropy.table import MaskedColumn, Table

from pipnick import log
from pipnick.datamodel.calibration import header_calibration_state
from pipnick.datamodel.dataset import DataSet
from pipnick.pkg.paths import calibration_file, step_dir
from pipnick.pkg.version import version


def astrometry_status(image, by_name):
    """
    Return the state of the astrometric solution of an image; see
    :mod:`pipnick.report.status`.

    A solution bootstrapped from another image is stale if the solution of
    that image has changed.

    Parameters
    ----------
    image : :class:`~pipnick.datamodel.image.ReducedImage`
        The image, which has been reduced.
    by_name : dict
        The science images of the data set, keyed by file name.

    Returns
    -------
    str
        ``'ok'``, ``'stale'``, or ``'missing'``.
    """
    hdr = image.header
    if 'PNASTMTH' not in hdr:
        return 'missing'
    ref = by_name.get(hdr.get('PNASTREF'))
    if ref is not None and ref.exists and ref.header.get('PNASTWH') != hdr.get('PNASTRH'):
        # Bootstrapped from an image whose solution changed
        return 'stale'
    return 'ok'


def photometry_status(image):
    """
    Return the state of the source table of an image; see
    :mod:`pipnick.report.status` and
    :meth:`~pipnick.datamodel.image.ReducedImage.photometry_state`.

    Parameters
    ----------
    image : :class:`~pipnick.datamodel.image.ReducedImage`
        The image, which has been reduced.

    Returns
    -------
    status : str
        ``'ok'``, ``'stale'`` (the image or its astrometric solution has
        changed), ``'missing'``, or ``'failed'`` (in the photometry summary).
    phot : :class:`astropy.table.Table`
        The source table, or None if it does not exist.
    """
    if not image.phot_file.is_file():
        failed = _summary_status(step_dir(image.dataset.rdxdir, 'photometry'), image.name)
        return 'failed' if failed == 'failed' else 'missing', None
    phot = Table.read(image.phot_file, format='ascii.ecsv')
    state, _ = image.photometry_state(phot)
    return 'ok' if state == 'ok' else 'stale', phot


def _summary_status(outdir, name):
    """
    Return the status of an image in the summary of a step, or None.
    """
    summary = outdir / 'summary.ecsv'
    if not summary.is_file():
        return None
    tbl = Table.read(summary, format='ascii.ecsv')
    if 'image' not in tbl.colnames:
        return None
    rows = tbl[tbl['image'] == name]
    return None if len(rows) == 0 else str(rows['status'][0])


def calibration_status(image):
    """
    Return the state of the photometric calibration in the header of an
    image; see :mod:`pipnick.report.status`.

    Parameters
    ----------
    image : :class:`~pipnick.datamodel.image.ReducedImage`
        The image, which has been reduced.

    Returns
    -------
    str
        ``'ok'``, ``'stale'`` (the calibration file has changed since it was
        applied), ``'missing'`` (not applied), or ``'n/a'`` (the night has no
        calibration file).
    """
    if image.night is None:
        return 'n/a'
    cfile = calibration_file(image.dataset.rdxdir, image.dataset.tag, image.night)
    if not cfile.is_file():
        return 'n/a'
    return header_calibration_state(image.header, cfile)


def mosaic_file(image, mosaic_dir=None):
    """
    Return the mosaic file of an image.

    Parameters
    ----------
    image : :class:`~pipnick.datamodel.image.ReducedImage`
        The image.
    mosaic_dir : str, :class:`pathlib.Path`, optional
        The directory with the mosaics.  If None, the ``mosaic`` subdirectory
        of the reduction directory is used.

    Returns
    -------
    :class:`pathlib.Path`
        The mosaic file, ``<mosaic>.fits``, or None if the image is not in a
        mosaic.
    """
    label = image.dataset.mosaic_labels(image.product)[0]
    if len(label) == 0:
        return None
    outdir = step_dir(image.dataset.rdxdir, 'mosaic') if mosaic_dir is None else Path(mosaic_dir)
    return outdir / f'{label}.fits'


def mosaic_status(image, mosaic_dir=None):
    """
    Return the state of the mosaic of an image; see
    :mod:`pipnick.report.status`.

    Only the image itself is checked, against its hashes in the ``INPUTS``
    table of the mosaic.  The mosaics also rebuild a mosaic whose settings,
    images, or WCS group have changed (``PNMOSHSH``; see
    :func:`~pipnick.pipelines.mosaic.mosaic_all`), which cannot be checked
    here, because a mosaic does not record all of its settings.

    Parameters
    ----------
    image : :class:`~pipnick.datamodel.image.ReducedImage`
        The image, which has been reduced.
    mosaic_dir : str, :class:`pathlib.Path`, optional
        The directory with the mosaics; see :func:`mosaic_file`.

    Returns
    -------
    str
        ``'ok'``, ``'stale'`` (the image or its solution has changed since the
        mosaic was constructed), ``'missing'``, or ``'n/a'`` (the image is not
        in a mosaic).
    """
    mfile = mosaic_file(image, mosaic_dir=mosaic_dir)
    if mfile is None:
        return 'n/a'
    if not mfile.is_file():
        return 'missing'
    with fits.open(mfile) as hdu:
        if 'INPUTS' not in hdu:
            return 'stale'
        inputs = Table(hdu['INPUTS'].data)
    rows = inputs[inputs['image'] == image.name]
    if (
        len(rows) == 0 or rows['image_hash'][0] != image.image_hash
        or rows['wcs_hash'][0] != image.solution_hash
    ):
        return 'stale'
    return 'ok'


def attention_flags(image, phot, camera):
    """
    Return the reasons an image needs attention, besides the state of its
    steps.

    Parameters
    ----------
    image : :class:`~pipnick.datamodel.image.ReducedImage`
        The image, which has been reduced.
    phot : :class:`astropy.table.Table`
        The source table of the image, or None.
    camera : :class:`~pipnick.cameras.cameras.NickelCamera`
        The camera, which sets the expected ranges; see
        :meth:`~pipnick.cameras.cameras.NickelCamera.expected_range`.

    Returns
    -------
    list
        The reasons: ``'version'`` (reduced by another pipnick version),
        ``'no bias'``, ``'no flat'``, ``'wcs rms'`` (an astrometric RMS above
        the expected range), and ``'fwhm'`` (a PSF FWHM outside it).
    """
    hdr = image.header
    flags = []
    if hdr.get('PNVERS') != version:
        flags.append('version')
    if hdr.get('PNBIAS', 'none') == 'none':
        flags.append('no bias')
    if hdr.get('PNFLAT', 'none') == 'none':
        flags.append('no flat')
    rms = camera.expected_range('astrometric_rms')
    if rms is not None and hdr.get('PNASTRMS') is not None and hdr['PNASTRMS'] > rms[1]:
        flags.append('wcs rms')
    fwhm = camera.expected_range('fwhm')
    value = None if phot is None else phot.meta.get('psf', {}).get('fwhm_arcsec')
    if fwhm is not None and value is not None and (value < fwhm[0] or value > fwhm[1]):
        flags.append('fwhm')
    return flags


def night_status(table, rdxdir=None, mosaic_dir=None):
    """
    Report the status of the processing of the science images of a data set.

    The status is determined from the products and their provenance (see
    :mod:`pipnick.report.status`); nothing is changed, and nothing is
    written.

    Parameters
    ----------
    table : str, :class:`pathlib.Path`
        A per-setup metadata file.
    rdxdir : str, :class:`pathlib.Path`, optional
        The reduction directory.  If None, the directory saved in the metadata
        file is used; see :func:`~pipnick.datamodel.dataset.reduction_dir`.
    mosaic_dir : str, :class:`pathlib.Path`, optional
        The directory with the mosaics; see :func:`mosaic_file`.

    Returns
    -------
    tbl : :class:`astropy.table.Table`
        One row per science image: its name, object, filter, exposure time
        and airmass; the state of the reduction (``reduced``), astrometric
        solution (``wcs``, with its RMS in arcsec), photometry (``phot``,
        with the number of sources and the PSF FWHM in arcsec), photometric
        calibration (``calibrated``, with the zero point) and mosaic
        (``mosaic``); and the reasons the image needs attention.  Unknown
        values are masked.
    next_step : str
        The next command to run.
    """
    dataset = DataSet.from_file(table, rdxdir=rdxdir, read_only=True)
    images = dataset.science_images(include_missing=True)
    by_name = {image.name: image for image in images}
    rows = []
    for image in images:
        values = image.values
        label = dataset.mosaic_labels(image.product)[0]
        row = {
            'image': image.name, 'object': values['object'], 'filter': values['filter'],
            'exptime': values['exptime'], 'airmass': values['airmass'], 'reduced': None,
            'wcs': 'n/a', 'rms': None, 'phot': 'n/a', 'nsrc': None, 'fwhm': None,
            'calibrated': 'n/a', 'zp': None, 'mosaic': 'n/a', 'attention': ''
        }
        rows.append(row)
        row['reduced'], _ = image.reduction_state()
        if row['reduced'] == 'missing':
            row['mosaic'] = 'n/a' if len(label) == 0 else 'missing'
            continue
        hdr = image.header
        row['wcs'] = astrometry_status(image, by_name)
        row['rms'] = hdr.get('PNASTRMS')
        row['phot'], phot = photometry_status(image)
        if phot is not None:
            row['nsrc'] = len(phot)
            row['fwhm'] = phot.meta.get('psf', {}).get('fwhm_arcsec')
        row['calibrated'] = calibration_status(image)
        row['zp'] = hdr.get('PNCALZP')
        row['mosaic'] = mosaic_status(image, mosaic_dir=mosaic_dir)
        row['attention'] = ', '.join(attention_flags(image, phot, dataset.camera))
    return _table(rows), _next_step(rows, dataset.metafile)


def _table(rows):
    """
    Construct the status table, masking the unknown values.
    """
    tbl = Table()
    columns = [
        'image', 'object', 'filter', 'exptime', 'airmass', 'reduced', 'wcs', 'rms', 'phot',
        'nsrc', 'fwhm', 'calibrated', 'zp', 'mosaic', 'attention'
    ]
    floats = {'exptime': '.1f', 'airmass': '.3f', 'rms': '.3f', 'fwhm': '.2f', 'zp': '.3f'}
    for name in columns:
        values = [row[name] for row in rows]
        missing = [v is None for v in values]
        if name in floats:
            tbl[name] = MaskedColumn(
                [0. if v is None else float(v) for v in values], mask=missing, dtype=float
            )
            tbl[name].info.format = floats[name]
        elif name == 'nsrc':
            tbl[name] = MaskedColumn(
                [0 if v is None else int(v) for v in values], mask=missing, dtype=int
            )
        else:
            tbl[name] = MaskedColumn(['' if v is None else str(v) for v in values], mask=missing)
    tbl.meta['pipnick_version'] = version
    return tbl


def _next_step(rows, metafile):
    """
    Return the next command to run: the first step with an image whose
    product is missing or stale.
    """
    steps = [
        ('reduced', 'pipnick_reduce'), ('wcs', 'pipnick_astrometry'),
        ('phot', 'pipnick_photometry'), ('calibrated', 'pipnick_calibrate'),
        ('mosaic', 'pipnick_mosaic')
    ]
    for column, script in steps:
        if any(row[column] in ['missing', 'stale', 'failed'] for row in rows):
            return f'{script} {metafile}'
    if all(row['calibrated'] == 'n/a' for row in rows):
        return (
            f'pipnick_calibrate {metafile}, if the night has standards; otherwise, every '
            'step is up to date.'
        )
    return 'nothing; every step is up to date.'


def log_status(tbl, next_step):
    """
    Log a status table and the next step; see :func:`night_status`.
    """
    log.info(
        'Status of the science images:\n'
        + '\n'.join(tbl.pformat(max_lines=-1, max_width=-1))
    )
    attention = tbl[np.asarray(tbl['attention'].filled('')) != '']
    if len(attention) > 0:
        names = list(attention['image'])
        listed = ', '.join(names[:5]) + (f', and {len(names) - 5} more' if len(names) > 5 else '')
        log.warning(f'{len(names)} image(s) need attention (see the attention column): {listed}')
    log.info(f'Next: {next_step}')
