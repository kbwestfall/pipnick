"""
Fit the photometric calibration of a night.

The calibration of each band is fit to the measurements of the standards in
the images of the night (see :mod:`pipnick.photometry.standards`) with
:func:`~pipnick.photometry.calibrate.fit_calibration`, with iterative outlier
rejection, a leave-one-out assessment of its accuracy, and a test for a drift
of the zero point; see :func:`calibrate_night`.
"""

import dataclasses
import datetime

import numpy as np

from astropy.table import MaskedColumn

from pipnick import log
from pipnick.datamodel.calibration import extinction_source
from pipnick.photometry.calibrate import fit_calibration, reject_outliers, residual_drift
from pipnick.photometry.landolt import catalog_path, landolt_band
from pipnick.photometry.standards import match_standards
from pipnick.pkg.exceptions import PipnickParameterError
from pipnick.pkg.provenance import file_hash
from pipnick.pkg.version import version


def _resolve_fixed(fix, band, terms, nominal):
    """
    Return the terms fixed by the user for one band, mapped to their values.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if a fixed term is not in the model, or the nominal value is
        requested for a term other than ``airmass`` or is not available.
    """
    fixed = {}
    for term, value in fix.items():
        if term not in terms:
            raise PipnickParameterError(f'Cannot fix {term}; it is not one of the terms.')
        if isinstance(value, dict):
            if band in value:
                fixed[term] = float(value[band])
        elif value == 'nominal':
            if term != 'airmass':
                raise PipnickParameterError('Only the airmass term has a nominal value.')
            if nominal is None:
                raise PipnickParameterError(
                    'The camera has no nominal extinction curve; the airmass term cannot be '
                    'fixed at its nominal value.'
                )
            fixed[term] = nominal[band]
        else:
            fixed[term] = float(value)
    return fixed



def calibrate_night(
    night, images, standards, catalog, camera, tag, settings, nominal, fix, input_hash
):
    """
    Fit the calibration of each band from the standards observed in a night.

    Parameters
    ----------
    night : str
        The night (UT date).
    images : list
        The images of the night; see :class:`~pipnick.datamodel.image.ReducedImage`.
    standards : list
        The images searched for standards.
    catalog : :class:`astropy.table.Table`
        The Landolt catalog.
    camera : :class:`~pipnick.cameras.cameras.NickelCamera`
        The camera.
    tag : str
        The setup tag.
    settings : :class:`~pipnick.pkg.settings.CalibrationSettings`
        The calibration settings.
    nominal : dict
        The nominal extinction of each band, or None.
    fix : dict
        The fixed terms; see :func:`~pipnick.pipelines.calibrate.parse_fix`.
    input_hash : str
        The hash of the inputs of the calibration, recorded in the
        calibration; see :func:`~pipnick.pipelines.calibrate.calibrate_all`.

    Returns
    -------
    content : dict
        The calibration; see :meth:`~pipnick.datamodel.calibration.Calibration.write`.
    meas : :class:`astropy.table.Table`
        The standard measurements and their residuals.
    """
    meas = match_standards(standards, catalog, settings)
    missing = np.ma.getmaskarray(meas['mag_cat_err']) & (meas['reason'] == '-')
    if np.any(missing):
        log.warning(
            f'{night}: {np.sum(missing)} usable standard measurements have no catalog error; '
            f'their errors include only the instrumental and systematic ({settings.sys_err} '
            'mag) errors.'
        )
    nmeas = len(meas)
    meas['mag_cal'] = MaskedColumn(np.zeros(nmeas), mask=np.ones(nmeas, dtype=bool))
    meas['resid'] = MaskedColumn(np.zeros(nmeas), mask=np.ones(nmeas, dtype=bool))
    meas['resid_err'] = MaskedColumn(np.zeros(nmeas), mask=np.ones(nmeas, dtype=bool))
    meas['resid_loo'] = MaskedColumn(np.zeros(nmeas), mask=np.ones(nmeas, dtype=bool))
    meas['used'] = np.zeros(nmeas, dtype=bool)
    meas['rejected'] = np.zeros(nmeas, dtype=bool)
    mjd0 = min(image.mjd for image in images)

    bands = sorted(
        {landolt_band(image.filter) for image in images if landolt_band(image.filter) is not None},
        key='UBVRI'.index
    )
    content = {
        'calibration': {
            'setup': tag, 'night': night, 'camera': camera.name,
            'catalog': 'Landolt (2009, 2013)', 'catalog_hash': file_hash(catalog_path()),
            'magnitude': settings.mag,
            'created': datetime.datetime.now(datetime.UTC).isoformat(timespec='seconds'),
            'pipnick_version': version, 'input_hash': input_hash
        },
        'band': {}
    }
    drift_sel = np.zeros(nmeas, dtype=bool)
    for band in bands:
        indx = np.where(meas['band'] == band)[0]
        bmeas = meas[indx]
        names = sorted({str(i) for i in bmeas['image']})
        fixed = _resolve_fixed(fix, band, settings.terms, nominal)
        nominal_k = None if nominal is None else nominal[band]
        if len(bmeas) == 0:
            result = {'status': 'failed', 'reason': 'no standards observed in this band'}
            content['band'][band] = _band_entry(
                band, settings.colors[band], bmeas, result, None, None, names
            )
            continue
        result = calibrate_band(bmeas, band, settings.terms, fixed, nominal_k, settings)
        for msg in result['warnings']:
            log.warning(f'{night} {msg}.')
        if result['status'] != 'ok':
            log.warning(f'{night}: the calibration of band {band} failed: {result["reason"]}.')
            content['band'][band] = _band_entry(
                band, settings.colors[band], bmeas, result, None, None, names
            )
            continue
        if result['fallback_reason'] is not None:
            log.warning(
                f'{night}: the extinction of band {band} was NOT MEASURED '
                f'({result["fallback_reason"]}); using the nominal value, '
                f'{nominal_k:.3f} mag/airmass.'
            )
        fit = result['fit']
        good = np.logical_not(
            np.ma.getmaskarray(bmeas['mag_inst']) | np.ma.getmaskarray(bmeas['mag_cat'])
            | np.ma.getmaskarray(bmeas['color'])
        )
        airmass = np.asarray(bmeas['airmass'], dtype=float)
        color = np.ma.filled(bmeas['color'], 0.)
        mag_cal = np.ma.filled(bmeas['mag_inst'], 0.) - fit.model(airmass, color)
        resid = mag_cal - np.ma.filled(bmeas['mag_cat'], 0.)
        resid_err = np.sqrt(np.asarray(bmeas['err'])**2 + fit.model_error(airmass, color)**2)
        meas['mag_cal'][indx] = np.ma.MaskedArray(mag_cal, mask=np.logical_not(good))
        meas['resid'][indx] = np.ma.MaskedArray(resid, mask=np.logical_not(good))
        meas['resid_err'][indx] = np.ma.MaskedArray(resid_err, mask=np.logical_not(good))
        loo = _leave_one_out(
            bmeas, band, settings.terms, fixed, nominal_k, settings, result['used']
        )
        meas['resid_loo'][indx] = loo
        meas['used'][indx] = result['used']
        meas['rejected'][indx] = result['rejected']
        drift = _drift(bmeas, result['used'], meas['resid'][indx], meas['resid_err'][indx], mjd0)
        _check_drift(f'{night} {band}', drift, settings.drift_warn)
        drift_sel[indx] = result['used']
        content['band'][band] = _band_entry(
            band, settings.colors[band], bmeas, result, loo, drift, names
        )

    drift = _drift(meas, drift_sel, meas['resid'], meas['resid_err'], mjd0)
    _check_drift(f'{night} (all bands)', drift, settings.drift_warn)
    if drift is not None:
        content['calibration']['drift'] = [drift[0], drift[1]]
    meas.meta = {
        'night': night, 'setup': tag, 'magnitude': settings.mag,
        'extinction': {
            band: extinction_source(entry) for band, entry in content['band'].items()
        },
        'residual': 'calibrated minus catalog magnitude'
    }
    for name in ['mag_cal', 'resid', 'resid_err', 'resid_loo']:
        meas[name].info.format = '.4f'
    return content, meas



def calibrate_band(meas, band, terms, fixed, nominal_k, settings, usable=None, strict=True):
    """
    Fit the calibration of one band.

    The fit is iterated: the measurements are checked (see below), the model
    is fit, and outliers are rejected (see
    :func:`~pipnick.photometry.calibrate.reject_outliers`), until no
    measurements are rejected or ``settings.maxiter`` rejection iterations
    are done.  Outliers are not rejected if this would leave no more
    measurements than fitted terms.

    Before each fit, the measurements are checked:

        - If there are fewer than ``settings.min_stars`` distinct stars,
          the calibration fails, unless ``settings.force`` is set, in
          which case a warning is issued.
        - If the airmass span is less than ``settings.min_airmass_span``
          while the ``airmass`` or ``airmass_color`` term is fit, the
          extinction cannot be measured.  If ``settings.fallback`` is set
          (the default), the extinction is fixed at its nominal value, the
          ``airmass_color`` term is fixed at 0, and the reason is recorded.
          Otherwise, the calibration fails.
        - If the color span is less than ``settings.min_color_span`` while
          a color term is fit, the calibration fails, unless
          ``settings.force`` is set.
        - If there are no more measurements than fitted terms, or the fitted
          terms are degenerate, the calibration fails.

    Parameters
    ----------
    meas : :class:`astropy.table.Table`
        The measurements of the band; see :func:`~pipnick.photometry.standards.match_standards`.
    band : str
        The band.
    terms : list
        The terms of the model.
    fixed : dict
        The terms fixed by the user, mapped to their values.
    nominal_k : float
        The nominal extinction coefficient of the band, or None if not
        available.
    settings : :class:`~pipnick.pkg.settings.CalibrationSettings`
        The calibration settings.
    usable : :class:`numpy.ndarray`, optional
        Boolean array selecting the measurements to use.  If None, the
        measurements without a reason not to use them are used.
    strict : bool, optional
        Apply the minimum-stars and color-span checks, and log the warnings.
        If False (used for the leave-one-out fits), only the fallback rule
        and the minimum number of measurements are applied.

    Returns
    -------
    result : dict
        ``status`` (``'ok'`` or ``'failed'``), ``reason`` (for a failure),
        ``fit`` (the :class:`~pipnick.photometry.calibrate.CalibrationFit`),
        ``source`` (the source of each term: ``'fit'``, ``'fixed'``, or
        ``'fallback'``), ``fallback_reason``, ``used`` and ``rejected``
        (boolean arrays), and ``warnings``.
    """
    _usable = np.asarray(meas['reason'] == '-') if usable is None else np.asarray(usable)
    used = _usable.copy()
    rejected = np.zeros(len(meas), dtype=bool)
    mag_c = np.ma.filled(meas['mag_cat'], 0.)
    color = np.ma.filled(meas['color'], 0.)
    mag_i = np.ma.filled(meas['mag_inst'], 0.)
    airmass = np.asarray(meas['airmass'], dtype=float)
    err = np.asarray(meas['err'], dtype=float)
    warnings = []
    result = {
        'status': 'failed', 'reason': None, 'fit': None, 'source': None,
        'fallback_reason': None, 'used': used, 'rejected': rejected, 'warnings': warnings
    }
    fit = None
    for iteration in range(settings.maxiter + 1):
        _fixed = dict(fixed)
        fallback = []
        nstars = len(np.unique(np.asarray(meas['star'])[used]))
        if strict and nstars < settings.min_stars:
            msg = f'only {nstars} standard stars (fewer than {settings.min_stars})'
            if not settings.force:
                result['reason'] = msg
                return result
            if f'{band}: {msg}' not in warnings:
                warnings.append(f'{band}: {msg}')
        span = float(np.ptp(airmass[used])) if np.any(used) else 0.
        if (
            any(t in terms and t not in _fixed for t in ['airmass', 'airmass_color'])
            and span < settings.min_airmass_span
        ):
            reason = f'airmass span {span:.2f} < {settings.min_airmass_span}'
            if not settings.fallback:
                result['reason'] = f'{reason}, and the extinction fallback is disabled'
                return result
            if 'airmass' in terms and 'airmass' not in _fixed:
                if nominal_k is None:
                    result['reason'] = f'{reason}, and there is no nominal extinction'
                    return result
                _fixed['airmass'] = nominal_k
                fallback.append('airmass')
            if 'airmass_color' in terms and 'airmass_color' not in _fixed:
                _fixed['airmass_color'] = 0.
                fallback.append('airmass_color')
            result['fallback_reason'] = reason
        cspan = float(np.ptp(color[used])) if np.any(used) else 0.
        if (
            strict and any(t in terms and t not in _fixed for t in ['color', 'color2'])
            and cspan < settings.min_color_span
        ):
            msg = f'color span {cspan:.2f} < {settings.min_color_span}'
            if not settings.force:
                result['reason'] = f'{msg}; fix the color terms (--fix)'
                return result
            if f'{band}: {msg}' not in warnings:
                warnings.append(f'{band}: {msg}')
        nfree = len([t for t in terms if t not in _fixed])
        nused = int(np.sum(used))
        if nused <= nfree:
            result['reason'] = f'only {nused} measurements for {nfree} fitted terms'
            return result
        try:
            fit = fit_calibration(
                mag_c, color, mag_i, airmass, terms, err=err, fixed=_fixed, gpm=used
            )
        except PipnickParameterError as e:
            result['reason'] = str(e)
            return result
        if iteration == settings.maxiter:
            break
        outliers = reject_outliers(fit, mag_c, color, mag_i, airmass, err=err, rej=settings.rej)
        if not np.any(outliers) or nused - int(np.sum(outliers)) <= nfree:
            break
        used = used & np.logical_not(outliers)
        rejected = rejected | outliers

    result.update(
        {
            'status': 'ok', 'fit': fit, 'used': used, 'rejected': rejected,
            'source': [
                'fallback' if t in fallback else ('fixed' if t in fixed else 'fit')
                for t in fit.terms
            ]
        }
    )
    if strict:
        for t, e, s, k in zip(fit.terms, fit.errors, result['source'], fit.coeff):
            if s == 'fit' and e > 0.1:
                warnings.append(f'{band}: the error of the {t} coefficient is large ({e:.3f})')
            if t == 'airmass' and s == 'fit' and k < 0:
                warnings.append(f'{band}: the fitted extinction is negative ({k:.3f})')
    return result



def _leave_one_out(meas, band, terms, fixed, nominal_k, settings, used):
    """
    Compute the leave-one-star-out residuals of a band.

    For each star used in the fit, the calibration is refit without any of
    its measurements, applying the same fallback rule (see
    :func:`calibrate_band`), and used to predict them.  The residuals (the
    calibrated minus the catalog magnitudes) are masked where the reduced fit
    has no remaining degrees of freedom or fails.
    """
    loo = np.ma.MaskedArray(np.zeros(len(meas)), mask=np.ones(len(meas), dtype=bool))
    stars = np.asarray(meas['star'])
    _settings = dataclasses.replace(settings, maxiter=0)
    for star in np.unique(stars[used]):
        others = used & (stars != star)
        result = calibrate_band(
            meas, band, terms, fixed, nominal_k, _settings, usable=others, strict=False
        )
        if result['status'] != 'ok' or result['fit'].dof == 0:
            continue
        indx = used & (stars == star)
        model = result['fit'].model(
            np.asarray(meas['airmass'])[indx], np.ma.filled(meas['color'], 0.)[indx]
        )
        loo[indx] = (
            np.ma.filled(meas['mag_inst'], 0.)[indx] - model
            - np.ma.filled(meas['mag_cat'], 0.)[indx]
        )
    return loo



def _drift(meas, sel, resid, resid_err, mjd0):
    """
    Return the drift of the residuals with time (mag/hour) and its error, or
    None; see :func:`~pipnick.photometry.calibrate.residual_drift`.
    """
    if np.sum(sel) == 0:
        return None
    hours = (np.asarray(meas['mjd'])[sel] - mjd0) * 24.
    return residual_drift(hours, np.ma.filled(resid[sel], 0.), np.ma.filled(resid_err[sel], 1.))



def _band_entry(band, color, meas, result, loo, drift, images):
    """
    Construct the entry of a band in the calibration file.
    """
    entry = {'status': result['status'], 'filter': band, 'color': color}
    if result['status'] != 'ok':
        entry['reason'] = result['reason']
        entry['n_stars'] = len(np.unique(np.asarray(meas['star'])[meas['reason'] == '-']))
        entry['n_meas'] = int(np.sum(meas['reason'] == '-'))
        entry['images'] = images
        return entry
    fit = result['fit']
    used = result['used']
    entry.update(
        {
            'terms': list(fit.terms), 'coefficients': [float(c) for c in fit.coeff],
            'errors': [float(e) for e in fit.errors], 'source': result['source']
        }
    )
    if result['fallback_reason'] is not None:
        entry['fallback_reason'] = result['fallback_reason']
    entry['covariance'] = [[float(c) for c in row] for row in fit.cov]
    entry['n_stars'] = len(np.unique(np.asarray(meas['star'])[used]))
    entry['n_meas'] = int(np.sum(used))
    entry['n_rejected'] = int(np.sum(result['rejected']))
    entry['n_unusable'] = int(np.sum(meas['reason'] != '-'))
    entry['rms'] = float(fit.rms)
    if fit.rchi2 is not None:
        entry['rchi2'] = float(fit.rchi2)
    n_loo = int(np.sum(np.logical_not(np.ma.getmaskarray(loo)) & used))
    entry['n_loo'] = len(
        np.unique(np.asarray(meas['star'])[np.logical_not(np.ma.getmaskarray(loo)) & used])
    )
    if n_loo > 0:
        entry['rms_loo'] = float(np.sqrt(np.mean(loo[used].compressed()**2)))
    if drift is not None:
        entry['drift'] = [drift[0], drift[1]]
    entry['airmass_range'] = [
        float(np.min(meas['airmass'][used])), float(np.max(meas['airmass'][used]))
    ]
    color_values = np.ma.filled(meas['color'], 0.)[used]
    entry['color_range'] = [float(np.min(color_values)), float(np.max(color_values))]
    entry['images'] = images
    return entry



def _check_drift(label, drift, limit):
    """
    Log the drift of the residuals, and warn if it is significant (more than
    3 sigma) and larger than ``limit`` (mag/hour).
    """
    if drift is None:
        return
    slope, err = drift
    msg = f'{label}: residual drift {slope:+.4f} +/- {err:.4f} mag/hour'
    if abs(slope) > 3 * err and abs(slope) > limit:
        log.warning(
            f'{msg}, which is significant; the zero point may have changed during the night.  '
            'The calibration does not model a drift.'
        )
    else:
        log.info(f'{msg}.')

