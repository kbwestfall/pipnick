r"""
Photometric calibration of the reduced science images with the Landolt
standard stars.

The Landolt standards (see :mod:`pipnick.photometry.landolt`) are found in the
photometry source tables by their celestial coordinates, so the astrometry
(:func:`~pipnick.pipelines.astrometry.astrometry_all`) and photometry
(:func:`~pipnick.pipelines.photometry.photometry_all`) must be run first.  The
calibration of each band is fit to the standards observed during one night
(see :mod:`pipnick.photometry.calibrate`) and saved to a TOML file that can be
edited.  Optionally, the calibration is applied to the source tables and
recorded in the image headers.

.. warning::

    If the standards of a band do not cover a sufficient range of airmass,
    the extinction cannot be measured, and its nominal value (see
    :mod:`pipnick.photometry.extinction`) is used instead.  This *fallback* is
    recorded in the calibration file, the image headers (``PNCALEXT``), the
    residual table, the QA figures, and the summary; see
    :func:`calibrate_all`.

.. note::

    The calibration does not include a drift of the zero point through the
    night.  The drift is measured from the residuals and reported, but not
    fit; fitting it is an area for improvement.
"""

import numpy as np

from astropy.table import Table

from pipnick import log
from pipnick.datamodel.calibration import Calibration
from pipnick.datamodel.dataset import DataSet
from pipnick.photometry.apply import apply_to_catalogs, apply_to_headers
from pipnick.photometry.calibrate import CALIBRATION_TERMS
from pipnick.photometry.extinction import nominal_extinction
from pipnick.photometry.landolt import catalog_path, landolt_band, load_landolt, parse_color
from pipnick.photometry.nightfit import calibrate_night
from pipnick.pkg.exceptions import (
    PipnickDataModelError, PipnickParameterError, PipnickPathError
)
from pipnick.pkg.paths import calibration_file, qa_figure, residuals_file, step_dir
from pipnick.pkg.provenance import file_hash, json_hash
from pipnick.pkg.report import StepReport
from pipnick.pkg.settings import CalibrationSettings
from pipnick.pkg.step import Attempt, needs_qa
from pipnick.pkg.summary import Summary
from pipnick.pkg.version import version
from pipnick.qa.calibration import plot_residuals


def _collect_images(dataset):
    """
    Collect the reduced science images with source tables.

    Images that do not exist, or have no source table or observation time,
    are skipped with a warning.  A warning is also issued for images whose
    filter cannot be calibrated (see
    :func:`~pipnick.photometry.landolt.landolt_band`), or whose source table
    has no celestial coordinates (see
    :attr:`~pipnick.datamodel.catalog.SourceCatalog.has_coordinates`).

    Parameters
    ----------
    dataset : :class:`~pipnick.datamodel.dataset.DataSet`
        The data set.

    Returns
    -------
    list
        The :class:`~pipnick.datamodel.image.ReducedImage` objects.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if an image is not from the setup, or its source table is out
        of date.
    """
    images = []
    unsupported = set()
    for image in dataset.science_images():
        if not image.phot_file.is_file():
            log.warning(f'{image.name} has no source table; run pipnick_photometry.  Skipping.')
            continue
        # Check that the source table is up to date
        image.photometry()
        if landolt_band(image.filter) is None and image.filter not in unsupported:
            unsupported.add(image.filter)
            log.warning(
                f'Filter {image.filter} cannot be calibrated: calibration currently requires '
                'UBVRI filters.'
            )
        if image.mjd is None:
            log.warning(f'{image.name} has no observation time; skipping.')
            continue
        if not image.photometry().has_coordinates:
            log.warning(
                f'{image.phot_file.name} has no celestial coordinates, so it cannot be searched '
                'for standards; run pipnick_astrometry, then pipnick_photometry.'
            )
        images.append(image)
    return images


def calibrate_all(
    table, rdxdir=None, images=None, terms=('zp', 'airmass', 'color'), fix=None, colors=None,
    mag='psf', match_radius=2., max_neighbor_flux=0.01, residual_radius=1.5, max_cat_err=0.05,
    min_cat_obs=1, sys_err=0.01, rej=3., maxiter=5, min_stars=3, min_airmass_span=0.2,
    min_color_span=0.3, force=False, fallback=True, drift_warn=0.02, qa=True, overwrite=False,
    apply=None, default_color=None, report=None
):
    r"""
    Calibrate the photometry of the reduced science images with the Landolt
    standards.

    The Landolt standards are found by matching the celestial coordinates of
    the sources in the photometry tables of all images to the catalog (see
    :func:`~pipnick.photometry.standards.match_standards`), so any image with
    a standard in its field is used, whatever its object.  The astrometry and
    photometry must have been run; images without celestial coordinates cannot
    be searched, and are skipped with a warning.

    Each band (filter U, B, V, R, or I; other filters cannot be calibrated)
    of each night is calibrated independently (see
    :func:`~pipnick.photometry.nightfit.calibrate_band`), with the model

    .. math::

        m - M = \sum_t k_t\ x_t;

    see :mod:`pipnick.photometry.calibrate`.  If the airmass span of the
    standards is too small to measure the extinction, the extinction is fixed
    at its nominal value (see :mod:`pipnick.photometry.extinction`), which is
    recorded in all of the products: the calibration file (the term's
    ``source`` is ``'fallback'``, with a ``fallback_reason`` and a comment),
    the image headers (``PNCALEXT``), the residual table, the QA figures, and
    the summary, which ends with a banner listing every fallback.

    The products, in the ``calibrate`` subdirectory of the reduction
    directory, are:

        - ``<setup>_<night>.toml``: the calibration, which may be edited; see
          :meth:`~pipnick.datamodel.calibration.Calibration.write`;
        - ``<setup>_<night>_residuals.ecsv``: the standard measurements, with
          their calibrated magnitudes, residuals (calibrated minus catalog),
          leave-one-star-out residuals (``resid_loo``; see
          :func:`~pipnick.photometry.nightfit.calibrate_night`), and whether
          they were used, rejected, or unusable (``reason``);
        - ``summary.ecsv``: the summary of the calibration of each band;
        - ``qa/<setup>_<night>_residuals.png``: the residuals against airmass,
          color, magnitude, and time (if ``qa`` is set and a band is
          calibrated).

    The residuals are also used to measure the drift of the zero point with
    time, which is reported but not fit.

    When rerun, a night is refit only if its inputs (the source tables and the
    settings) have changed, or ``overwrite`` is set.  A calibration file that
    has been edited is kept (and used to apply the calibration) unless
    ``overwrite`` is set, with a warning if its inputs have changed.

    The errors are statistical only: fixed and fallback terms have no error,
    and the systematic error they introduce is not included.

    Parameters
    ----------
    table : str, :class:`pathlib.Path`
        The per-setup metadata file used to reduce the images; see
        :func:`~pipnick.pipelines.reduction.reduce_all`.
    rdxdir : str, :class:`pathlib.Path`, optional
        The reduction directory.  If None, the directory saved in the metadata
        file is used.
    images : list, optional
        Only search these images for standards, given by file name, file name
        without the extension, or combination group (``comb_id``).  This does
        not limit the images the calibration is applied to.
    terms : list, optional
        The terms of the model; see
        :data:`~pipnick.photometry.calibrate.CALIBRATION_TERMS`.
    fix : dict, optional
        Terms to fix; see :func:`parse_fix`.
    colors : dict, optional
        The color used to calibrate each band, e.g., ``{'B': 'B-V'}``; the
        defaults are given by
        :func:`~pipnick.photometry.landolt.default_colors`.
    mag : str, optional
        The instrumental magnitudes to calibrate: ``'psf'`` (``mag_inst``) or
        ``'aper'`` (from ``flux_aper``).  No aperture correction is applied,
        so a calibration of aperture magnitudes is valid only for magnitudes
        measured in the same aperture.
    match_radius : float, optional
        The radius (arcsec) for matching sources to the standards.
    max_neighbor_flux : float, optional
        Standards with a neighbor within their aperture brighter than this
        fraction of their flux are not used; see
        :meth:`~pipnick.datamodel.catalog.SourceCatalog.neighbors`.
    residual_radius : float, optional
        Sources detected in the residuals of a standard (in a later iteration
        of the PSF photometry) within this radius (in FWHM) do not exclude it;
        see :meth:`~pipnick.datamodel.catalog.SourceCatalog.neighbors`.
    max_cat_err : float, optional
        Standards with a catalog error (in the band) larger than this (mag)
        are not used.
    min_cat_obs : int, optional
        Standards observed fewer times than this by Landolt are not used.
    sys_err : float, optional
        A systematic error (mag) added in quadrature to the error of each
        measurement.
    rej : float, optional
        The outlier rejection threshold, in standard deviations.
    maxiter : int, optional
        The maximum number of rejection iterations.
    min_stars : int, optional
        The minimum number of standard stars in each band.
    min_airmass_span : float, optional
        The minimum airmass span required to fit the extinction.
    min_color_span : float, optional
        The minimum color span (mag) required to fit the color terms.
    force : bool, optional
        Downgrade the minimum-stars and color-span errors to warnings.
    fallback : bool, optional
        If the airmass span is too small, fix the extinction at its nominal
        value.  If False, the calibration of the band fails instead.
    drift_warn : float, optional
        Warn if the residuals drift by more than this (mag/hour), and the
        drift is significant.
    qa : bool, optional
        Write the QA figures: a figure is redrawn when a night is fit, and
        drawn from its residual table for a night that is not refit but lacks
        one.
    overwrite : bool, optional
        Refit, even if the inputs have not changed or the calibration file
        has been edited.
    apply : str, optional
        Apply the calibration: ``'catalog'`` (see
        :func:`~pipnick.photometry.apply.apply_to_catalogs`), ``'image'``
        (see :func:`~pipnick.photometry.apply.apply_to_headers`), or
        ``'both'``.  It is
        applied to all of the science images of each night whose band is
        calibrated.
    default_color : dict, optional
        The color assumed for sources without a measured color, keyed by the
        color, e.g., ``{'B-V': 0.6}``.  Colors not given default to 0.
    report : :class:`~pipnick.pkg.report.StepReport`, optional
        Collects the problems of the calibration for its closing verdict: the
        bands that could not be calibrated, and those whose extinction was not
        measured.  If None, a new one is used.

    Returns
    -------
    list
        The paths to the calibration files.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickPathError`
        Raised if the reduction directory does not exist.
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if no source table has celestial coordinates, a source table
        is out of date, or no band is calibrated.
    """
    _report = StepReport('CALIBRATION') if report is None else report
    settings = CalibrationSettings(
        terms=terms, colors=colors, mag=mag, match_radius=match_radius,
        max_neighbor_flux=max_neighbor_flux, residual_radius=residual_radius,
        max_cat_err=max_cat_err, min_cat_obs=min_cat_obs, sys_err=sys_err, rej=rej,
        maxiter=maxiter, min_stars=min_stars, min_airmass_span=min_airmass_span,
        min_color_span=min_color_span, force=force, fallback=fallback, drift_warn=drift_warn
    )
    if apply not in [None, 'catalog', 'image', 'both']:
        raise PipnickParameterError(f'Unknown apply option {apply}.')
    _fix = {} if fix is None else fix
    _default_color = {} if default_color is None else default_color
    for color in _default_color:
        parse_color(color)
    dataset = DataSet.from_file(table, rdxdir=rdxdir)
    _report.check_version(dataset.metafile.name, dataset.table.meta.get('pipnick_version'))
    _rdxdir = dataset.rdxdir
    tag = dataset.tag
    camera = dataset.camera
    if not _rdxdir.is_dir():
        raise PipnickPathError(f'{_rdxdir} does not exist; reduce the data first.')
    if settings.mag == 'aper':
        log.warning(
            'Calibrating the aperture magnitudes, without an aperture correction: the '
            'calibration is valid only for magnitudes measured in the same aperture.'
        )
    cal_dir = step_dir(_rdxdir, 'calibrate')
    log.info(f'Photometric calibration of the images in {_rdxdir}')
    for band, color in settings.colors.items():
        if band not in parse_color(color):
            log.warning(f'The color {color} used for band {band} does not include it.')
    nominal = nominal_extinction(camera, ['U', 'B', 'V', 'R', 'I'])
    if nominal is None:
        log.warning(f'Camera {camera.name} has no nominal extinction curve.')
    else:
        log.info(
            'Nominal extinction (mag/airmass): '
            + ', '.join(f'{b} {k:.3f}' for b, k in nominal.items())
        )

    selected = {
        p.comb_id
        for p in (dataset.products('science') if images is None else dataset.select(images))
    }
    all_images = _collect_images(dataset)
    if not any(image.photometry().has_coordinates for image in all_images):
        raise PipnickDataModelError(
            'No source table has celestial coordinates, so no standards can be found; run '
            'pipnick_astrometry, then pipnick_photometry.'
        )
    catalog = load_landolt()

    outputs = []
    summary = []
    nights = sorted({image.night for image in all_images})
    for night in nights:
        failed = {'night': night, 'band': '-', 'status': 'failed', 'action': 'failed'}
        attempt = Attempt(
            _report, f'The calibration of {night}', row=failed, then='the other nights'
        )
        with attempt:
            night_images = [image for image in all_images if image.night == night]
            standards = [
                image for image in night_images
                if image.photometry().has_coordinates and landolt_band(image.filter) is not None
                and image.airmass is not None and image.product.comb_id in selected
            ]
            input_hash = json_hash(
                {
                    'settings': settings.as_dict(),
                    'fix': _fix, 'nominal': nominal, 'catalog': file_hash(catalog_path()),
                    'tables': [
                        [
                            image.phot_file.name, image.photometry().meta.get('settings_hash'),
                            image.photometry().meta.get('image_hash'),
                            image.photometry().meta.get('wcs_hash')
                        ]
                        for image in standards
                    ],
                    'bands': sorted({str(landolt_band(i.filter)) for i in night_images})
                }
            )
            ofile = calibration_file(_rdxdir, tag, night)
            rfile = residuals_file(_rdxdir, tag, night)
            qa_file = qa_figure(cal_dir, rfile.stem)
            content = None
            if ofile.is_file() and not overwrite:
                calibration = Calibration.read(ofile)
                existing, edited = calibration.content, calibration.edited
                same_inputs = existing['calibration'].get('input_hash') == input_hash
                if edited:
                    log.warning(
                        f'{ofile.name} has been edited; keeping it (use --overwrite to refit).'
                    )
                    if not same_inputs:
                        log.warning(
                            f'The inputs of {ofile.name} have changed since it was fit; the '
                            'edited calibration may hide changes in the photometry or '
                            'settings.'
                        )
                    content = existing
                    action = 'edited (kept)'
                    _report.note(f'{ofile.name} has been edited, and was kept.')
                elif same_inputs:
                    log.info(f'{ofile.name} is up to date.')
                    content = existing
                    action = 'up to date'
                if content is not None:
                    _report.check_version(
                        ofile.name, content['calibration'].get('pipnick_version')
                    )
            if content is None:
                if len(standards) == 0:
                    log.warning(
                        f'{night}: no images with celestial coordinates (and an airmass) to '
                        'search for standards.'
                    )
                content, meas = calibrate_night(
                    night, night_images, standards, catalog, camera, tag, settings, nominal, _fix,
                    input_hash
                )
                Calibration(content).write(ofile)
                log.info(f'Wrote calibrate/{ofile.name}')
                meas.meta['pipnick_version'] = version
                meas.write(rfile, format='ascii.ecsv', overwrite=True)
                log.info(f'Wrote calibrate/{rfile.name}')
                if qa:
                    _write_qa(meas, content, qa_file, _report)
                action = 'fit'
            elif needs_qa(qa, False, [qa_file]) and rfile.is_file():
                _write_qa(Table.read(rfile), content, qa_file, _report)
            outputs.append(ofile)
            summary += Calibration(content).summary_rows(night, action)

            if apply in ['image', 'both']:
                apply_to_headers(night_images, content, ofile, report=_report)
            if apply in ['catalog', 'both']:
                apply_to_catalogs(night_images, content, ofile, settings, _default_color)
        if attempt.failed:
            summary.append(failed)

    _summarize(summary, cal_dir)
    for row in summary:
        if row['status'] != 'ok':
            _report.attention(
                f'{row["night"]} {row["band"]}: the calibration failed: {row["reason"]}',
                warn=False
            )
        elif row['extinction'] == 'fallback':
            _report.attention(
                f'{row["night"]} {row["band"]}: the extinction was not measured; the nominal '
                'value was used (see the banner above).', warn=False
            )
    nok = sum(row['status'] == 'ok' for row in summary)
    text = f'{nok} of {len(summary)} bands calibrated; summary in calibrate/summary.ecsv.'
    if nok == 0:
        _report.fail('No band could be calibrated; see the summary above.')
        _report.finish(text)
        raise PipnickDataModelError('No band could be calibrated; see the summary above.')
    following = None
    if 'mosaic' in dataset.table.colnames and np.any(dataset.table['mosaic'] != ''):
        following = f'pipnick_mosaic {dataset.metafile}'
    _report.finish(text, next_step=following)
    return outputs


def _write_qa(meas, content, ofile, report):
    """
    Write the QA figure of the calibration of a night; see
    :func:`~pipnick.qa.calibration.plot_residuals`.  A failure is recorded in
    the report, and does not stop the calibration.
    """
    attempt = Attempt(report, f'The QA figure {ofile.name}', problem=True)
    with attempt:
        plot_residuals(meas, content, ofile)
    if attempt.failed:
        return
    if ofile.is_file():
        log.info(f'Wrote calibrate/qa/{ofile.name}')


def _summarize(summary, cal_dir):
    """
    Log the summary of the calibration of each band, write it to
    ``calibrate/summary.ecsv``, and end with a banner listing the bands whose
    extinction was not measured.

    The summary gives, for each night and band, the status of the
    calibration, the numbers of stars, measurements, and rejected
    measurements, the rms and leave-one-out rms of the residuals, how the
    extinction was found and its value, the zero point above the atmosphere
    and its error, the drift of the residuals, what was done with the
    calibration file (``action``: ``fit``, ``up to date``, or ``edited
    (kept)``), and the reason for a failure or fallback; see
    :meth:`~pipnick.datamodel.calibration.Calibration.summary_rows`.
    """
    columns = [
        'night', 'band', 'status', 'n_stars', 'n_meas', 'n_rejected', 'rms', 'rms_loo',
        'extinction', 'k_airmass', 'zp', 'zp_err', 'drift', 'drift_err', 'action', 'reason'
    ]
    floats = ['rms', 'rms_loo', 'k_airmass', 'zp', 'zp_err', 'drift', 'drift_err']
    table = Summary(
        'Calibration', columns,
        dtypes={'n_stars': int, 'n_meas': int, 'n_rejected': int} | {k: float for k in floats},
        formats={k: '.3f' for k in floats},
        logged=[c for c in columns if c not in ['zp_err', 'drift_err', 'reason']]
    )
    for row in summary:
        table.add(**row)
    tbl = table.write(cal_dir)
    if tbl is None:
        log.warning('No bands to calibrate.')
        return
    for row in tbl:
        if row['status'] != 'ok':
            log.warning(f'{row["night"]} {row["band"]}: calibration failed: {row["reason"]}')
    fallback = tbl[tbl['extinction'] == 'fallback']
    if len(fallback) > 0:
        lines = [
            '*' * 78,
            'EXTINCTION NOT MEASURED: the nominal extinction was used for these bands, so',
            'their calibration does not include a measurement of the extinction on the night,',
            'and their errors exclude its systematic error:'
        ]
        lines += [f'    {r["night"]} {r["band"]}: {r["reason"]}' for r in fallback]
        lines += ['*' * 78]
        log.warning('\n' + '\n'.join(lines))


def parse_assignments(values, kind):
    """
    Parse command-line assignments of the form ``key=value``.

    Parameters
    ----------
    values : list
        The assignments.
    kind : str
        What is assigned, used in error messages.

    Returns
    -------
    dict
        The values (strings), keyed by the keys.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if an assignment is not of the form ``key=value``.
    """
    out = {}
    for item in [] if values is None else values:
        key, sep, value = str(item).partition('=')
        if sep != '=' or len(key) == 0 or len(value) == 0:
            raise PipnickParameterError(f'Invalid {kind} {item}; must be of the form key=value.')
        out[key.strip()] = value.strip()
    return out



def parse_fix(values):
    """
    Parse the ``--fix`` command-line option.

    Parameters
    ----------
    values : list
        Each entry is ``<term>=<value>``, where the value is ``nominal`` (only
        for ``airmass``; see :mod:`pipnick.photometry.extinction`), a number
        used for all bands, or a comma-separated list of ``<band>:<number>``
        pairs, e.g. ``airmass=B:0.28,V:0.17``.

    Returns
    -------
    dict
        The value of each fixed term: ``'nominal'``, a float, or a dictionary
        with the float for each band.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if an entry is invalid.
    """
    fix = {}
    for term, value in parse_assignments(values, 'fixed term').items():
        if term not in CALIBRATION_TERMS:
            raise PipnickParameterError(
                f'Cannot fix unknown term {term}; options are {", ".join(CALIBRATION_TERMS)}.'
            )
        if value == 'nominal':
            fix[term] = value
            continue
        try:
            if ':' in value:
                fix[term] = {
                    b.split(':')[0].strip(): float(b.split(':')[1]) for b in value.split(',')
                }
            else:
                fix[term] = float(value)
        except (ValueError, IndexError) as e:
            raise PipnickParameterError(f'Invalid value for fixed term {term}: {value}.') from e
    return fix

