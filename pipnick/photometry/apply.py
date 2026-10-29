"""
Apply the photometric calibration of a night to its reduced images and
source tables.
"""

import numpy as np

from astropy.io import fits
from astropy.table import MaskedColumn, Table

from pipnick import log
from pipnick.datamodel.calibration import Calibration, header_calibration_state
from pipnick.photometry.calibrate import iterate_colors, zero_point
from pipnick.photometry.landolt import landolt_band, parse_color
from pipnick.pkg.paths import calibrated_table
from pipnick.pkg.provenance import file_hash
from pipnick.pkg.version import version
from pipnick.utils.crossmatch import link_sources


def apply_to_headers(images, content, cfile, report=None):
    r"""
    Record the calibration in the headers of the reduced images.

    The pixels are not changed.  For each image whose band is calibrated, the
    keywords give the calibrated magnitude of a source with count rate
    :math:`r` (electrons per second) and color :math:`C`:

    .. math::

        M = -2.5\log_{10}(r) + {\rm PNCALZP} + {\rm PNCALCT}\ C
            + {\rm PNCALC2}\ C^2;

    see :func:`~pipnick.photometry.calibrate.zero_point`.  The keywords are:

        - ``PNCALZP``, ``PNCALZE``: the zero point at the airmass of the image
          and its statistical error;
        - ``PNCALCT``, ``PNCALC2``: the color coefficients;
        - ``PNCALCOL``: the color;
        - ``PNCALMAG``: the instrumental magnitude calibrated (``psf`` or
          ``aper``);
        - ``PNCALSYS``: the magnitude system (``Vega``);
        - ``PNCALEXT``: how the extinction was set: ``fit``, ``fixed``,
          ``fallback`` (the nominal value; not measured), or ``none``;
        - ``PNCALFIL``, ``PNCALHSH``: the calibration file and its hash;
        - ``PHOTZP``: the zero point for :math:`C = 0`, equal to ``PNCALZP``.

    The keywords are written only if ``PNCALHSH`` differs from the hash of the
    calibration file.  The keywords are removed from images whose band is not
    calibrated.

    Parameters
    ----------
    images : list
        The images; see :class:`~pipnick.datamodel.image.ReducedImage`.
    content : dict
        The calibration; see :meth:`~pipnick.datamodel.calibration.Calibration.write`.
    cfile : :class:`pathlib.Path`
        The calibration file.
    report : :class:`~pipnick.pkg.report.StepReport`, optional
        Records the version (``PNCALVER``) of the calibrations in the headers
        that are up to date.

    Returns
    -------
    list
        The names of the images that could not be calibrated.
    """
    calibrations = Calibration(content).bands()
    chash = file_hash(cfile)
    skipped = []
    for image in images:
        path = image.path
        cal = calibrations.get(landolt_band(image.filter))
        if cal is None or image.airmass is None:
            skipped.append(path.name)
            log.warning(
                f'{path.name}: no calibration for filter {image.photometry().meta.get("filter")}'
                + ('' if cal is None else ' (no airmass)') + '; its header is not calibrated.'
            )
            stale = [
                k for k in fits.getheader(path) if k.startswith('PNCAL') or k == 'PHOTZP'
            ]
            if len(stale) > 0:
                with fits.open(path, mode='update') as hdu:
                    for key in stale:
                        del hdu[0].header[key]
                log.warning(f'{path.name}: removed the out-of-date calibration keywords.')
            continue
        hdr = fits.getheader(path)
        if header_calibration_state(hdr, cfile) == 'ok':
            log.info(f'{path.name}: header calibration is up to date.')
            if report is not None:
                report.check_version(f'{path.name} (calibration)', hdr.get('PNCALVER'))
            continue
        zp = zero_point(cal['terms'], cal['coeff'], cal['cov'], float(image.airmass))
        cards = [
            ('PNCALZP', zp['zp'], 'Zero point (e-/s) at this airmass, color 0'),
            ('PNCALZE', zp['zp_err'], 'Statistical error of PNCALZP'),
            ('PNCALCT', zp['c1'], 'Coefficient of the color'),
            ('PNCALC2', zp['c2'], 'Coefficient of the squared color'),
            ('PNCALCOL', cal['color_name'], 'Color used by the calibration'),
            ('PNCALMAG', content['calibration']['magnitude'], 'Instrumental magnitude type'),
            ('PNCALSYS', 'Vega', 'Magnitude system (Landolt)'),
            ('PNCALEXT', cal['extinction'], 'Extinction: fit/fixed/fallback/none'),
            ('PNCALFIL', cfile.name, 'Calibration file'),
            ('PNCALHSH', chash, 'Hash of the calibration file'),
            ('PNCALVER', version, 'pipnick version (calibration)'),
            ('PHOTZP', zp['zp'], 'Zero point (e-/s) for color 0')
        ]
        with fits.open(path, mode='update') as hdu:
            for key, value, comment in cards:
                hdu[0].header[key] = (value, comment)
        log.info(
            f'{path.name}: recorded the calibration in the header (zero point '
            f'{zp["zp"]:.3f}, extinction {cal["extinction"]}).'
        )
    return skipped



def apply_to_catalogs(images, content, cfile, settings, default_color):
    """
    Calibrate the source tables of the images of a night.

    The sources in the images of each object are identified by their
    celestial coordinates (see :func:`~pipnick.utils.crossmatch.link_sources`),
    and calibrated with :func:`~pipnick.photometry.calibrate.iterate_colors`,
    so that the color of each source is computed from its calibrated
    magnitudes, averaged over the images in each band.  Sources without a
    measured color use the default color, and are flagged
    (``color_assumed``).

    The calibrated table of each image is written to
    ``calibrate/<image>_cal.ecsv``; the source tables of the photometry are
    not changed.  Its columns are the source ``id``, ``x``, ``y``, ``ra``,
    ``dec``, the instrumental magnitude calibrated (``mag_inst`` and
    ``mag_inst_err``), the calibrated magnitude (``mag``) and its statistical
    error (``mag_err``), ``color_used``, ``color_assumed``, and the identifier
    of the source across the images of the object (``source_id``).

    Parameters
    ----------
    images : list
        The images; see :class:`~pipnick.datamodel.image.ReducedImage`.
    content : dict
        The calibration; see :meth:`~pipnick.datamodel.calibration.Calibration.write`.
    cfile : :class:`pathlib.Path`
        The calibration file.
    settings : :class:`~pipnick.pkg.settings.CalibrationSettings`
        The calibration settings.
    default_color : dict
        The default colors, keyed by the color, e.g., ``{'B-V': 0.6}``.

    Returns
    -------
    list
        The paths to the calibrated tables.
    """
    calibrations = Calibration(content).bands()
    defaults = {parse_color(c): float(v) for c, v in default_color.items()}
    mag_type = content['calibration']['magnitude']
    chash = file_hash(cfile)
    outputs = []
    usable = []
    for image in images:
        if landolt_band(image.filter) not in calibrations or image.airmass is None:
            log.warning(
                f'{image.name}: its band is not calibrated; not calibrating its source '
                'table.'
            )
            continue
        usable.append(image)
    objects = sorted({str(image.photometry().meta.get('object')) for image in usable})
    for obj in objects:
        group = [image for image in usable if str(image.photometry().meta.get('object')) == obj]
        # Images without celestial coordinates are not matched
        ra = [
            np.ma.asarray(i.photometry()['ra'], dtype=float) if i.photometry().has_coordinates
            else np.ma.masked_all(len(i.photometry())) for i in group
        ]
        dec = [
            np.ma.asarray(i.photometry()['dec'], dtype=float) if i.photometry().has_coordinates
            else np.ma.masked_all(len(i.photometry())) for i in group
        ]
        ids = link_sources(ra, dec, settings.match_radius)
        mags = []
        errs = []
        for image in group:
            m, e = image.photometry().instrumental_mag(mag_type)
            mags.append(m)
            errs.append(e)
        n = [len(m) for m in mags]
        mag, mag_err, color, assumed = iterate_colors(
            np.ma.concatenate(mags), np.ma.concatenate(errs),
            np.concatenate([np.full(k, float(i.airmass)) for k, i in zip(n, group)]),
            np.concatenate(ids),
            np.concatenate([[landolt_band(i.filter)] * k for k, i in zip(n, group)]),
            calibrations, default_color=defaults
        )
        start = 0
        for image, k, sid, m, e in zip(group, n, ids, mags, errs):
            phot = image.photometry()
            cal = calibrations[landolt_band(image.filter)]
            tbl = Table()
            for name in ['id', 'x', 'y', 'ra', 'dec']:
                tbl[name] = phot[name]
            tbl['mag_inst'] = MaskedColumn(m, unit='mag')
            tbl['mag_inst_err'] = MaskedColumn(e, unit='mag')
            tbl['mag'] = MaskedColumn(mag[start:start + k], unit='mag')
            tbl['mag_err'] = MaskedColumn(mag_err[start:start + k], unit='mag')
            tbl['color_used'] = MaskedColumn(color[start:start + k], unit='mag')
            tbl['color_assumed'] = assumed[start:start + k]
            tbl['source_id'] = sid
            for name in ['mag_inst', 'mag_inst_err', 'mag', 'mag_err', 'color_used']:
                tbl[name].info.format = '.4f'
            tbl.meta = {
                'image': image.name, 'object': obj, 'band': landolt_band(image.filter),
                'color': cal['color_name'], 'magnitude': mag_type, 'system': 'Vega',
                'extinction': cal['extinction'], 'calibration': cfile.name,
                'calibration_hash': chash, 'airmass': image.airmass,
                'errors': 'statistical only; fixed and fallback terms have no error',
                'pipnick_version': version
            }
            ofile = calibrated_table(image.path)
            tbl.write(ofile, format='ascii.ecsv', overwrite=True)
            outputs.append(ofile)
            log.info(
                f'Wrote {ofile.name}: {np.sum(np.logical_not(np.ma.getmaskarray(tbl["mag"])))} '
                f'calibrated sources, {np.sum(tbl["color_assumed"])} with an assumed color.'
            )
            start += k
    return outputs

