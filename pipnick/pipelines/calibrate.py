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

import datetime
import hashlib
import json
import tomllib
from pathlib import Path

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import MaskedColumn, Table
from astropy.time import Time
from matplotlib.figure import Figure

from pipnick import log
from pipnick.photometry.calibrate import (
    CALIBRATION_TERMS, check_terms, fit_calibration, iterate_colors, reject_outliers,
    residual_drift, zero_point
)
from pipnick.photometry.extinction import nominal_extinction
from pipnick.photometry.landolt import (
    catalog_color, catalog_path, load_landolt, match_landolt, parse_color
)
from pipnick.pipelines.photometry import image_hash, select_products, solution_hash
from pipnick.pipelines.reduction import load_setup, plan_products
from pipnick.pkg.exceptions import (
    PipnickDataModelError, PipnickParameterError, PipnickPathError
)


def default_colors():
    """
    Return the default color used to calibrate each band: U-B for U, B-V for
    B and V, V-R for R, and V-I for I.
    """
    return {'U': 'U-B', 'B': 'B-V', 'V': 'B-V', 'R': 'V-R', 'I': 'V-I'}


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


def _hash(content):
    """
    Return the first 16 hexadecimal digits of the SHA-1 hash of the JSON
    representation of ``content``.
    """
    return hashlib.sha1(json.dumps(content, sort_keys=True).encode()).hexdigest()[:16]


def _file_hash(path):
    """
    Return the first 16 hexadecimal digits of the SHA-1 hash of a file.
    """
    return hashlib.sha1(Path(path).read_bytes()).hexdigest()[:16]


def _night(mjd):
    """
    Return the UT date (``YYYY-MM-DD``) of an MJD.
    """
    return Time(mjd, format='mjd', scale='utc').to_value('iso', subfmt='date')


def _band(filtername):
    """
    Return the Landolt band of a filter, or None if it is not one of UBVRI.
    """
    return filtername if filtername in ['U', 'B', 'V', 'R', 'I'] else None


def _instrumental(phot, mag):
    r"""
    Return the instrumental magnitudes of a source table and their errors.

    For ``mag='psf'``, these are ``mag_inst`` and ``mag_inst_err``.  For
    ``mag='aper'``, they are :math:`-2.5\log_{10}(f / t)` for the aperture
    flux :math:`f` and the exposure time :math:`t`, masked where the flux is
    masked or not positive.
    """
    if mag == 'psf':
        return (
            np.ma.asarray(phot['mag_inst'], dtype=float),
            np.ma.asarray(phot['mag_inst_err'], dtype=float)
        )
    flux = np.ma.asarray(phot['flux_aper'], dtype=float)
    ferr = np.ma.asarray(phot['flux_aper_err'], dtype=float)
    bad = np.ma.getmaskarray(flux) | (flux.filled(0.) <= 0)
    _flux = np.where(bad, 1., flux.filled(1.))
    value = -2.5 * np.log10(_flux / float(phot.meta['exptime']))
    err = 2.5 / np.log(10.) * ferr.filled(0.) / _flux
    return (
        np.ma.MaskedArray(value, mask=bad),
        np.ma.MaskedArray(err, mask=bad | np.ma.getmaskarray(ferr))
    )


def neighbors(phot, index, residual_radius=1.5, max_neighbor_flux=0.01):
    """
    Classify the sources near a source in a photometry table.

    The PSF photometry of bright stars often detects faint, spurious sources
    in the residuals of their fits, very close to the star, which are then fit
    together with it (so its ``group_size`` is larger than 1).  These
    *residual detections* are identified as sources detected in a later
    iteration than the source (``iter_detected``) and within
    ``residual_radius`` FWHM of it; they are counted, but do not exclude the
    source.  Any other source within the photometric aperture of the source
    (``aper_radius_pix`` in the table metadata) with more than
    ``max_neighbor_flux`` of its PSF flux is a *contaminating* neighbor.

    Parameters
    ----------
    phot : :class:`astropy.table.Table`
        The photometry source table; see
        :func:`~pipnick.pipelines.photometry.photometry_all`.
    index : int
        The row of the source.
    residual_radius : float, optional
        The radius, in FWHM, within which later detections are taken to be in
        the residuals of the source.
    max_neighbor_flux : float, optional
        The fraction of the PSF flux of the source above which a neighbor
        within the aperture contaminates it.

    Returns
    -------
    n_residual, n_contaminating : int
        The numbers of residual detections and contaminating neighbors.
    """
    x = np.asarray(phot['x'], dtype=float)
    y = np.asarray(phot['y'], dtype=float)
    dist = np.hypot(x - x[index], y - y[index])
    flux = np.ma.filled(np.ma.asarray(phot['flux_psf'], dtype=float), 0.)
    iteration = (
        np.asarray(phot['iter_detected']) if 'iter_detected' in phot.colnames
        else np.ones(len(phot), dtype=int)
    )
    other = np.arange(len(phot)) != index
    residual = (
        other & (iteration > iteration[index])
        & (dist < residual_radius * phot.meta['psf']['fwhm_pix'])
    )
    contaminating = (
        other & np.logical_not(residual) & (dist < phot.meta['aper_radius_pix'])
        & (flux > max_neighbor_flux * flux[index])
    )
    return int(np.sum(residual)), int(np.sum(contaminating))


def match_standards(images, catalog, settings):
    """
    Find the Landolt standards in the source tables of a night.

    Parameters
    ----------
    images : list
        The images, as returned by :func:`_collect_images`.
    catalog : :class:`astropy.table.Table`
        The Landolt catalog; see :func:`~pipnick.photometry.landolt.load_landolt`.
    settings : dict
        The calibration settings, including the color used to calibrate each
        band (``colors``); see :func:`calibrate_all`.

    Returns
    -------
    :class:`astropy.table.Table`
        One row per matched standard and image, with the image, band, star,
        its number of Landolt observations, the match separation (arcsec), the
        airmass and MJD, the catalog magnitude and color and their errors, the
        instrumental magnitude and its error, the number of residual
        detections near the standard (``n_residual``; see :func:`neighbors`),
        the measurement error used to weight the fit (``err``), and the reason
        the measurement is not used (``'-'`` if it is usable).
    """
    rows = []
    for image in images:
        phot = image['phot']
        ra = np.ma.asarray(phot['ra'], dtype=float)
        dec = np.ma.asarray(phot['dec'], dtype=float)
        good = np.logical_not(np.ma.getmaskarray(ra) | np.ma.getmaskarray(dec))
        epoch = Time(image['mjd'], format='mjd', scale='utc').jyear
        isrc, icat, sep = match_landolt(
            catalog, ra[good].filled(0.), dec[good].filled(0.), epoch,
            radius=settings['match_radius']
        )
        isrc = np.where(good)[0][isrc]
        band = image['band']
        mag_i, err_i = _instrumental(phot, settings['mag'])
        cat = catalog[icat]
        mag_c = np.ma.asarray(cat[band], dtype=float)
        err_c = np.ma.asarray(cat[f'e_{band}'], dtype=float)
        color, color_err = catalog_color(cat, settings['colors'][band])
        flags = np.asarray(phot['flags'])
        aper_masked = np.ma.getmaskarray(np.ma.asarray(phot['flux_aper']))
        for j, (s, c) in enumerate(zip(isrc, icat)):
            n_residual, n_contaminating = neighbors(
                phot, s, residual_radius=settings['residual_radius'],
                max_neighbor_flux=settings['max_neighbor_flux']
            )
            if np.ma.getmaskarray(mag_c)[j]:
                reason = 'no catalog magnitude'
            elif np.ma.getmaskarray(color)[j]:
                reason = 'no catalog color'
            elif np.ma.getmaskarray(mag_i)[s]:
                reason = 'masked magnitude'
            elif flags[s] != 0:
                reason = 'flagged'
            elif aper_masked[s]:
                reason = 'saturated or otherwise masked'
            elif n_contaminating > 0:
                reason = 'contaminated'
            elif (
                not np.ma.getmaskarray(err_c)[j] and err_c[j] > settings['max_cat_err']
            ):
                reason = 'catalog error'
            elif cat['nobs'][j] < settings['min_cat_obs']:
                reason = 'few catalog observations'
            else:
                reason = '-'
            rows.append(
                {
                    'image': image['path'].stem, 'band': band, 'star': str(cat['name'][j]),
                    'nobs': int(cat['nobs'][j]), 'source_id': int(phot['id'][s]),
                    'sep': float(sep[j]), 'airmass': float(image['airmass']),
                    'mjd': float(image['mjd']), 'mag_cat': mag_c[j], 'mag_cat_err': err_c[j],
                    'color': color[j], 'color_err': color_err[j], 'mag_inst': mag_i[s],
                    'mag_inst_err': err_i[s], 'n_residual': n_residual, 'reason': reason
                }
            )
    return _measurement_table(rows, settings['sys_err'])


def _measurement_table(rows, sys_err):
    """
    Construct the table of standard measurements; see
    :func:`match_standards`.
    """
    names = {
        'image': str, 'band': str, 'star': str, 'nobs': int, 'source_id': int, 'sep': float,
        'airmass': float, 'mjd': float, 'mag_cat': float, 'mag_cat_err': float,
        'color': float, 'color_err': float, 'mag_inst': float, 'mag_inst_err': float,
        'n_residual': int, 'reason': str
    }
    tbl = Table()
    for name, dtype in names.items():
        values = [r[name] for r in rows]
        if dtype is float:
            mask = [v is np.ma.masked for v in values]
            data = np.array([0. if m else float(v) for v, m in zip(values, mask)], dtype=float)
            tbl[name] = MaskedColumn(data, mask=np.array(mask, dtype=bool))
        else:
            tbl[name] = (
                np.array(values, dtype=dtype) if len(values) > 0
                else np.array([], dtype=dtype if dtype is int else 'U1')
            )
    # The error of each measurement of m - M: the instrumental and catalog
    # errors (0 if unknown) and the systematic error, in quadrature
    tbl['err'] = np.sqrt(
        np.ma.filled(tbl['mag_inst_err'], 0.)**2 + np.ma.filled(tbl['mag_cat_err'], 0.)**2
        + sys_err**2
    )
    formatted = [
        'sep', 'airmass', 'mag_cat', 'mag_cat_err', 'color', 'color_err', 'mag_inst',
        'mag_inst_err', 'err'
    ]
    for name in formatted:
        tbl[name].info.format = '.4f'
    tbl['sep'].unit = 'arcsec'
    return tbl


def calibrate_band(meas, band, terms, fixed, nominal_k, settings, usable=None, strict=True):
    """
    Fit the calibration of one band.

    The fit is iterated: the measurements are checked (see below), the model
    is fit, and outliers are rejected (see
    :func:`~pipnick.photometry.calibrate.reject_outliers`), until no
    measurements are rejected or ``settings['maxiter']`` rejection iterations
    are done.  Outliers are not rejected if this would leave no more
    measurements than fitted terms.

    Before each fit, the measurements are checked:

        - If there are fewer than ``settings['min_stars']`` distinct stars,
          the calibration fails, unless ``settings['force']`` is set, in
          which case a warning is issued.
        - If the airmass span is less than ``settings['min_airmass_span']``
          while the ``airmass`` or ``airmass_color`` term is fit, the
          extinction cannot be measured.  If ``settings['fallback']`` is set
          (the default), the extinction is fixed at its nominal value, the
          ``airmass_color`` term is fixed at 0, and the reason is recorded.
          Otherwise, the calibration fails.
        - If the color span is less than ``settings['min_color_span']`` while
          a color term is fit, the calibration fails, unless
          ``settings['force']`` is set.
        - If there are no more measurements than fitted terms, or the fitted
          terms are degenerate, the calibration fails.

    Parameters
    ----------
    meas : :class:`astropy.table.Table`
        The measurements of the band; see :func:`match_standards`.
    band : str
        The band.
    terms : list
        The terms of the model.
    fixed : dict
        The terms fixed by the user, mapped to their values.
    nominal_k : float
        The nominal extinction coefficient of the band, or None if not
        available.
    settings : dict
        The calibration settings; see :func:`calibrate_all`.
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
    for iteration in range(settings['maxiter'] + 1):
        _fixed = dict(fixed)
        fallback = []
        nstars = len(np.unique(np.asarray(meas['star'])[used]))
        if strict and nstars < settings['min_stars']:
            msg = f'only {nstars} standard stars (fewer than {settings["min_stars"]})'
            if not settings['force']:
                result['reason'] = msg
                return result
            if f'{band}: {msg}' not in warnings:
                warnings.append(f'{band}: {msg}')
        span = float(np.ptp(airmass[used])) if np.any(used) else 0.
        if (
            any(t in terms and t not in _fixed for t in ['airmass', 'airmass_color'])
            and span < settings['min_airmass_span']
        ):
            reason = f'airmass span {span:.2f} < {settings["min_airmass_span"]}'
            if not settings['fallback']:
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
            and cspan < settings['min_color_span']
        ):
            msg = f'color span {cspan:.2f} < {settings["min_color_span"]}'
            if not settings['force']:
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
        if iteration == settings['maxiter']:
            break
        outliers = reject_outliers(fit, mag_c, color, mag_i, airmass, err=err, rej=settings['rej'])
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
    _settings = dict(settings, maxiter=0)
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


def _toml_value(value):
    """
    Format a value for the calibration file.
    """
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return repr(float(value))
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, (list, tuple)):
        return '[' + ', '.join(_toml_value(v) for v in value) + ']'
    raise PipnickParameterError(f'Cannot write {type(value)} to the calibration file.')


def _toml_text(content):
    """
    Construct the text of a calibration file; see :func:`write_calibration`.
    """
    cal = content['calibration']
    lines = [
        '# Photometric calibration written by pipnick_calibrate.  The coefficients (and the',
        '# terms) may be edited; an edited file is kept, and used by --apply, unless',
        '# pipnick_calibrate is run with --overwrite.  The model of each band is',
        '#     m_inst - M = sum_t k_t x_t,',
        '# for the terms x_t (zp: 1; airmass: X; color: C; color2: C^2; airmass_color: X C).',
        '# The errors are statistical only: fixed and fallback terms have no error, and',
        '# the systematic error they introduce is not included.'
    ]
    if cal['magnitude'] == 'aper':
        lines += [
            '# WARNING: calibrated with aperture magnitudes, without an aperture correction;',
            '# the calibration is valid only for magnitudes measured in the same aperture.'
        ]
    lines += ['', '[calibration]']
    lines += [f'{k} = {_toml_value(v)}' for k, v in cal.items() if k != 'fit_hash']
    if 'fit_hash' in cal:
        lines.append(f'fit_hash = {_toml_value(cal["fit_hash"])}')
    for band, entry in content.get('band', {}).items():
        lines.append('')
        if entry['status'] != 'ok':
            lines.append(f'# CALIBRATION FAILED for this band: {entry["reason"]}')
        elif 'fallback' in entry['source']:
            lines += [
                '# WARNING: extinction NOT MEASURED for this band; the nominal value was used',
                f'# ({entry["fallback_reason"]}).'
            ]
        lines.append(f'[band.{band}]')
        lines += [f'{k} = {_toml_value(v)}' for k, v in entry.items()]
    return '\n'.join(lines) + '\n'


def write_calibration(path, content):
    """
    Write a calibration file.

    The file is TOML, with a ``[calibration]`` table (the setup, night,
    catalog, magnitude, creation time, and the hash of the inputs) and a
    ``[band.<band>]`` table for each band; see :func:`calibrate_all`.  The
    ``fit_hash`` entry is the hash of the content as written (see
    :func:`calibration_hash`), used to detect edits.

    Parameters
    ----------
    path : :class:`pathlib.Path`
        The output file.
    content : dict
        The calibration, without ``fit_hash``.
    """
    content['calibration'].pop('fit_hash', None)
    # The hash is computed from the content as it is read back, so that it
    # does not depend on how the values are formatted
    parsed = tomllib.loads(_toml_text(content))
    content['calibration']['fit_hash'] = calibration_hash(parsed)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_toml_text(content))


def calibration_hash(content):
    """
    Return the hash of the content of a calibration file, excluding its
    ``fit_hash`` entry.
    """
    _content = json.loads(json.dumps(content))
    _content.get('calibration', {}).pop('fit_hash', None)
    return _hash(_content)


def read_calibration(path):
    """
    Read a calibration file.

    Returns
    -------
    content : dict
        The calibration.
    edited : bool
        True if the content no longer matches its ``fit_hash``, i.e., it has
        been edited.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the file cannot be parsed, or the coefficients of a band do
        not match its terms.
    """
    try:
        content = tomllib.loads(Path(path).read_text())
    except tomllib.TOMLDecodeError as e:
        raise PipnickDataModelError(f'Cannot parse {path}: {e}') from e
    for band, entry in content.get('band', {}).items():
        if entry.get('status') != 'ok':
            continue
        entry['terms'] = check_terms(entry['terms'])
        n = len(entry['terms'])
        if len(entry['coefficients']) != n:
            raise PipnickDataModelError(
                f'{path}: band {band} has {len(entry["coefficients"])} coefficients for {n} '
                'terms.'
            )
        cov = np.asarray(entry.get('covariance', np.zeros((n, n))), dtype=float)
        if cov.shape != (n, n):
            raise PipnickDataModelError(f'{path}: the covariance of band {band} is not {n}x{n}.')
    return content, calibration_hash(content) != content['calibration'].get('fit_hash')


def _collect_images(products, tag, phot_dir):
    """
    Collect the reduced science images and their source tables.

    Images that do not exist, or have no source table, are skipped with a
    warning.  The other images are returned as dictionaries with the path,
    the source table, the band (None if the filter is not one of UBVRI), the
    MJD, airmass, and night, and whether the table has celestial coordinates.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if an image is not from the setup, or its source table is out
        of date.
    """
    images = []
    unsupported = set()
    for product in products:
        path = product['file']
        pfile = phot_dir / f'{path.stem}_phot.ecsv'
        if not path.is_file():
            log.warning(f'{path.name} does not exist; reduce the data first.  Skipping.')
            continue
        if not pfile.is_file():
            log.warning(f'{path.name} has no source table; run pipnick_photometry.  Skipping.')
            continue
        hdr = fits.getheader(path)
        if hdr.get('PNSETUP') != tag:
            raise PipnickDataModelError(f'{path.name} is not from setup {tag}.')
        phot = Table.read(pfile, format='ascii.ecsv')
        if phot.meta.get('image_hash') != image_hash(hdr):
            raise PipnickDataModelError(
                f'The source table {pfile.name} is out of date: {path.name} has changed since '
                'it was measured.  Rerun pipnick_photometry.'
            )
        if phot.meta.get('wcs_hash') != solution_hash(hdr):
            raise PipnickDataModelError(
                f'The celestial coordinates in {pfile.name} are out of date: the astrometric '
                f'solution of {path.name} has changed.  Rerun pipnick_photometry to update '
                'them.'
            )
        filtername = phot.meta.get('filter')
        band = None if filtername is None else _band(str(filtername))
        if band is None and filtername not in unsupported:
            unsupported.add(filtername)
            log.warning(
                f'Filter {filtername} cannot be calibrated: calibration currently requires '
                'UBVRI filters.'
            )
        mjd = phot.meta.get('mjd')
        if mjd is None:
            log.warning(f'{path.name} has no observation time; skipping.')
            continue
        coords = (
            phot.meta.get('wcs_method', 'none') != 'none'
            and np.any(np.logical_not(np.ma.getmaskarray(np.ma.asarray(phot['ra']))))
        )
        if not coords:
            log.warning(
                f'{pfile.name} has no celestial coordinates, so it cannot be searched for '
                'standards; run pipnick_astrometry, then pipnick_photometry.'
            )
        images.append(
            {
                'path': path, 'phot_file': pfile, 'phot': phot, 'band': band, 'mjd': float(mjd),
                'airmass': phot.meta.get('airmass'), 'night': _night(float(mjd)),
                'coords': coords, 'product': product
            }
        )
    return images


def calibrate_night(night, images, standards, catalog, camera, tag, settings, nominal, fix):
    """
    Fit the calibration of each band from the standards observed in a night.

    Parameters
    ----------
    night : str
        The night (UT date).
    images : list
        The images of the night; see :func:`_collect_images`.
    standards : list
        The images searched for standards.
    catalog : :class:`astropy.table.Table`
        The Landolt catalog.
    camera : :class:`~pipnick.cameras.cameras.NickelCamera`
        The camera.
    tag : str
        The setup tag.
    settings : dict
        The calibration settings; see :func:`calibrate_all`.
    nominal : dict
        The nominal extinction of each band, or None.
    fix : dict
        The fixed terms; see :func:`parse_fix`.

    Returns
    -------
    content : dict
        The calibration; see :func:`write_calibration`.
    meas : :class:`astropy.table.Table`
        The standard measurements and their residuals.
    """
    meas = match_standards(standards, catalog, settings)
    missing = np.ma.getmaskarray(meas['mag_cat_err']) & (meas['reason'] == '-')
    if np.any(missing):
        log.warning(
            f'{night}: {np.sum(missing)} usable standard measurements have no catalog error; '
            f'their errors include only the instrumental and systematic ({settings["sys_err"]} '
            'mag) errors.'
        )
    nmeas = len(meas)
    meas['mag_cal'] = MaskedColumn(np.zeros(nmeas), mask=np.ones(nmeas, dtype=bool))
    meas['resid'] = MaskedColumn(np.zeros(nmeas), mask=np.ones(nmeas, dtype=bool))
    meas['resid_err'] = MaskedColumn(np.zeros(nmeas), mask=np.ones(nmeas, dtype=bool))
    meas['resid_loo'] = MaskedColumn(np.zeros(nmeas), mask=np.ones(nmeas, dtype=bool))
    meas['used'] = np.zeros(nmeas, dtype=bool)
    meas['rejected'] = np.zeros(nmeas, dtype=bool)
    mjd0 = min(image['mjd'] for image in images)

    bands = sorted(
        {image['band'] for image in images if image['band'] is not None},
        key='UBVRI'.index
    )
    content = {
        'calibration': {
            'setup': tag, 'night': night, 'camera': camera.name,
            'catalog': 'Landolt (2009, 2013)', 'catalog_hash': _file_hash(catalog_path()),
            'magnitude': settings['mag'],
            'created': datetime.datetime.now(datetime.UTC).isoformat(timespec='seconds'),
            'input_hash': settings['input_hash']
        },
        'band': {}
    }
    drift_sel = np.zeros(nmeas, dtype=bool)
    for band in bands:
        indx = np.where(meas['band'] == band)[0]
        bmeas = meas[indx]
        names = sorted({str(i) for i in bmeas['image']})
        fixed = _resolve_fixed(fix, band, settings['terms'], nominal)
        nominal_k = None if nominal is None else nominal[band]
        if len(bmeas) == 0:
            result = {'status': 'failed', 'reason': 'no standards observed in this band'}
            content['band'][band] = _band_entry(
                band, settings['colors'][band], bmeas, result, None, None, names
            )
            continue
        result = calibrate_band(bmeas, band, settings['terms'], fixed, nominal_k, settings)
        for msg in result['warnings']:
            log.warning(f'{night} {msg}.')
        if result['status'] != 'ok':
            log.warning(f'{night}: the calibration of band {band} failed: {result["reason"]}.')
            content['band'][band] = _band_entry(
                band, settings['colors'][band], bmeas, result, None, None, names
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
            bmeas, band, settings['terms'], fixed, nominal_k, settings, result['used']
        )
        meas['resid_loo'][indx] = loo
        meas['used'][indx] = result['used']
        meas['rejected'][indx] = result['rejected']
        drift = _drift(bmeas, result['used'], meas['resid'][indx], meas['resid_err'][indx], mjd0)
        _check_drift(f'{night} {band}', drift, settings['drift_warn'])
        drift_sel[indx] = result['used']
        content['band'][band] = _band_entry(
            band, settings['colors'][band], bmeas, result, loo, drift, names
        )

    drift = _drift(meas, drift_sel, meas['resid'], meas['resid_err'], mjd0)
    _check_drift(f'{night} (all bands)', drift, settings['drift_warn'])
    if drift is not None:
        content['calibration']['drift'] = [drift[0], drift[1]]
    meas.meta = {
        'night': night, 'setup': tag, 'magnitude': settings['mag'],
        'extinction': {
            band: _extinction_source(entry) for band, entry in content['band'].items()
        },
        'residual': 'calibrated minus catalog magnitude'
    }
    for name in ['mag_cal', 'resid', 'resid_err', 'resid_loo']:
        meas[name].info.format = '.4f'
    return content, meas


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


def _extinction_source(entry):
    """
    Return how the extinction of a band was set: ``'fit'``, ``'fixed'``,
    ``'fallback'``, ``'none'`` (no airmass term), or ``'failed'``.
    """
    if entry.get('status') != 'ok':
        return 'failed'
    if 'airmass' not in entry['terms']:
        return 'none'
    return entry['source'][entry['terms'].index('airmass')]


def calibrate_all(
    table, rdxdir=None, images=None, terms=('zp', 'airmass', 'color'), fix=None, colors=None,
    mag='psf', match_radius=2., max_neighbor_flux=0.01, residual_radius=1.5, max_cat_err=0.05,
    min_cat_obs=1, sys_err=0.01, rej=3., maxiter=5, min_stars=3, min_airmass_span=0.2,
    min_color_span=0.3, force=False, fallback=True, drift_warn=0.02, qa=False, overwrite=False,
    apply=None, default_color=None
):
    r"""
    Calibrate the photometry of the reduced science images with the Landolt
    standards.

    The Landolt standards are found by matching the celestial coordinates of
    the sources in the photometry tables of all images to the catalog (see
    :func:`match_standards`), so any image with a standard in its field is
    used, whatever its object.  The astrometry and photometry must have been
    run; images without celestial coordinates cannot be searched, and are
    skipped with a warning.

    Each band (filter U, B, V, R, or I; other filters cannot be calibrated)
    of each night is calibrated independently (see :func:`calibrate_band`),
    with the model

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
          :func:`write_calibration`;
        - ``<setup>_<night>_residuals.ecsv``: the standard measurements, with
          their calibrated magnitudes, residuals (calibrated minus catalog),
          leave-one-star-out residuals (``resid_loo``; see
          :func:`_leave_one_out`), and whether they were used, rejected, or
          unusable (``reason``);
        - ``summary.ecsv``: the summary of the calibration of each band;
        - ``qa/<setup>_<night>_residuals.png``: the residuals against airmass,
          color, magnitude, and time (if ``qa`` is set).

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
        defaults are given by :func:`default_colors`.
    mag : str, optional
        The instrumental magnitudes to calibrate: ``'psf'`` (``mag_inst``) or
        ``'aper'`` (from ``flux_aper``).  No aperture correction is applied,
        so a calibration of aperture magnitudes is valid only for magnitudes
        measured in the same aperture.
    match_radius : float, optional
        The radius (arcsec) for matching sources to the standards.
    max_neighbor_flux : float, optional
        Standards with a neighbor within their aperture brighter than this
        fraction of their flux are not used; see :func:`neighbors`.
    residual_radius : float, optional
        Sources detected in the residuals of a standard (in a later iteration
        of the PSF photometry) within this radius (in FWHM) do not exclude it;
        see :func:`neighbors`.
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
        Write the QA figures.
    overwrite : bool, optional
        Refit, even if the inputs have not changed or the calibration file
        has been edited.
    apply : str, optional
        Apply the calibration: ``'catalog'`` (see :func:`apply_to_catalogs`),
        ``'image'`` (see :func:`apply_to_headers`), or ``'both'``.  It is
        applied to all of the science images of each night whose band is
        calibrated.
    default_color : dict, optional
        The color assumed for sources without a measured color, keyed by the
        color, e.g., ``{'B-V': 0.6}``.  Colors not given default to 0.

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
    metafile, metadata, _rdxdir, tag, camera = load_setup(table, rdxdir=rdxdir)
    if not _rdxdir.is_dir():
        raise PipnickPathError(f'{_rdxdir} does not exist; reduce the data first.')
    if mag not in ['psf', 'aper']:
        raise PipnickParameterError(f'Unknown magnitude {mag}; must be psf or aper.')
    if apply not in [None, 'catalog', 'image', 'both']:
        raise PipnickParameterError(f'Unknown apply option {apply}.')
    if mag == 'aper':
        log.warning(
            'Calibrating the aperture magnitudes, without an aperture correction: the '
            'calibration is valid only for magnitudes measured in the same aperture.'
        )
    cal_dir = _rdxdir / 'calibrate'
    log.info(f'Photometric calibration of the images in {_rdxdir}')

    _colors = default_colors()
    _colors.update({} if colors is None else colors)
    for band, color in _colors.items():
        bands = parse_color(color)
        if band not in bands:
            log.warning(f'The color {color} used for band {band} does not include it.')
    _fix = {} if fix is None else fix
    _default_color = {} if default_color is None else default_color
    for color in _default_color:
        parse_color(color)
    settings = {
        'terms': check_terms(list(terms)), 'colors': _colors, 'mag': mag,
        'match_radius': float(match_radius), 'max_neighbor_flux': float(max_neighbor_flux),
        'residual_radius': float(residual_radius), 'max_cat_err': float(max_cat_err),
        'min_cat_obs': int(min_cat_obs), 'sys_err': float(sys_err), 'rej': float(rej),
        'maxiter': int(maxiter), 'min_stars': int(min_stars),
        'min_airmass_span': float(min_airmass_span), 'min_color_span': float(min_color_span),
        'force': bool(force), 'fallback': bool(fallback), 'drift_warn': float(drift_warn)
    }
    nominal = nominal_extinction(camera, ['U', 'B', 'V', 'R', 'I'])
    if nominal is None:
        log.warning(f'Camera {camera.name} has no nominal extinction curve.')
    else:
        log.info(
            'Nominal extinction (mag/airmass): '
            + ', '.join(f'{b} {k:.3f}' for b, k in nominal.items())
        )

    products = [p for p in plan_products(metadata, _rdxdir) if p['frametype'] == 'science']
    selected = {
        p['comb_id'] for p in (products if images is None else select_products(products, images))
    }
    all_images = _collect_images(products, tag, _rdxdir / 'photometry')
    if not any(image['coords'] for image in all_images):
        raise PipnickDataModelError(
            'No source table has celestial coordinates, so no standards can be found; run '
            'pipnick_astrometry, then pipnick_photometry.'
        )
    catalog = load_landolt()

    outputs = []
    summary = []
    nights = sorted({image['night'] for image in all_images})
    for night in nights:
        night_images = [image for image in all_images if image['night'] == night]
        standards = [
            image for image in night_images
            if image['coords'] and image['band'] is not None and image['airmass'] is not None
            and image['product']['comb_id'] in selected
        ]
        settings['input_hash'] = _hash(
            {
                'settings': {k: v for k, v in settings.items() if k != 'input_hash'},
                'fix': _fix, 'nominal': nominal, 'catalog': _file_hash(catalog_path()),
                'tables': [
                    [
                        image['phot_file'].name, image['phot'].meta.get('settings_hash'),
                        image['phot'].meta.get('image_hash'), image['phot'].meta.get('wcs_hash')
                    ]
                    for image in standards
                ],
                'bands': sorted({str(i['band']) for i in night_images})
            }
        )
        ofile = cal_dir / f'{tag}_{night}.toml'
        content = None
        if ofile.is_file() and not overwrite:
            existing, edited = read_calibration(ofile)
            same_inputs = existing['calibration'].get('input_hash') == settings['input_hash']
            if edited:
                log.warning(
                    f'{ofile.name} has been edited; keeping it (use --overwrite to refit).'
                )
                if not same_inputs:
                    log.warning(
                        f'The inputs of {ofile.name} have changed since it was fit; the edited '
                        'calibration may hide changes in the photometry or settings.'
                    )
                content = existing
                status = 'edited (kept)'
            elif same_inputs:
                log.info(f'{ofile.name} is up to date.')
                content = existing
                status = 'up to date'
        if content is None:
            if len(standards) == 0:
                log.warning(
                    f'{night}: no images with celestial coordinates (and an airmass) to search '
                    'for standards.'
                )
            content, meas = calibrate_night(
                night, night_images, standards, catalog, camera, tag, settings, nominal, _fix
            )
            write_calibration(ofile, content)
            log.info(f'Wrote {ofile}')
            rfile = cal_dir / f'{tag}_{night}_residuals.ecsv'
            meas.write(rfile, format='ascii.ecsv', overwrite=True)
            log.info(f'Wrote {rfile}')
            if qa:
                _write_qa(meas, content, cal_dir / 'qa' / f'{tag}_{night}_residuals.png')
            status = 'fit'
        outputs.append(ofile)
        summary += _summary_rows(night, content, status)

        if apply in ['image', 'both']:
            apply_to_headers(night_images, content, ofile)
        if apply in ['catalog', 'both']:
            apply_to_catalogs(night_images, content, ofile, cal_dir, settings, _default_color)

    _summarize(summary, cal_dir)
    if not any(row['status'] == 'ok' for row in summary):
        raise PipnickDataModelError('No band could be calibrated; see the summary above.')
    return outputs


def _calibrations(content):
    """
    Return the calibrations of the bands of a calibration file, in the form
    used by :func:`~pipnick.photometry.calibrate.iterate_colors`.
    """
    calibrations = {}
    for band, entry in content.get('band', {}).items():
        if entry.get('status') != 'ok':
            continue
        n = len(entry['terms'])
        calibrations[band] = {
            'terms': entry['terms'], 'coeff': np.asarray(entry['coefficients'], dtype=float),
            'cov': np.asarray(entry.get('covariance', np.zeros((n, n))), dtype=float),
            'color': parse_color(entry['color']), 'color_name': entry['color'],
            'extinction': _extinction_source(entry)
        }
    return calibrations


def apply_to_headers(images, content, cfile):
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
        The images; see :func:`_collect_images`.
    content : dict
        The calibration; see :func:`write_calibration`.
    cfile : :class:`pathlib.Path`
        The calibration file.

    Returns
    -------
    list
        The names of the images that could not be calibrated.
    """
    calibrations = _calibrations(content)
    chash = _file_hash(cfile)
    skipped = []
    for image in images:
        path = image['path']
        cal = calibrations.get(image['band'])
        if cal is None or image['airmass'] is None:
            skipped.append(path.name)
            log.warning(
                f'{path.name}: no calibration for filter {image["phot"].meta.get("filter")}'
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
        if hdr.get('PNCALHSH') == chash:
            log.info(f'{path.name}: header calibration is up to date.')
            continue
        zp = zero_point(cal['terms'], cal['coeff'], cal['cov'], float(image['airmass']))
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


def _cross_match(images, radius):
    """
    Identify the same sources in images of one object, by their celestial
    coordinates.

    The sources of each image (in order) are matched to the sources found so
    far, within ``radius`` (arcsec); unmatched sources are added.  Sources in
    images without celestial coordinates are not matched.

    Returns
    -------
    list
        The source identifier of each source in each image.
    """
    ref_id = []
    ref_ra = []
    ref_dec = []
    ids = []
    next_id = 0
    for image in images:
        phot = image['phot']
        sid = np.arange(next_id, next_id + len(phot))
        next_id += len(phot)
        ids.append(sid)
        if not image['coords'] or len(phot) == 0:
            continue
        ra = np.ma.asarray(phot['ra'], dtype=float)
        dec = np.ma.asarray(phot['dec'], dtype=float)
        good = np.where(np.logical_not(np.ma.getmaskarray(ra) | np.ma.getmaskarray(dec)))[0]
        new = good
        if len(ref_id) > 0 and len(good) > 0:
            ref = SkyCoord(ref_ra, ref_dec, unit='deg')
            coo = SkyCoord(ra[good].filled(0.), dec[good].filled(0.), unit='deg')
            near, sep, _ = coo.match_to_catalog_sky(ref)
            matched = sep.arcsec < radius
            sid[good[matched]] = np.asarray(ref_id)[near[matched]]
            new = good[np.logical_not(matched)]
        ref_id += [int(s) for s in sid[new]]
        ref_ra += [float(r) for r in ra[new]]
        ref_dec += [float(d) for d in dec[new]]
    return ids


def apply_to_catalogs(images, content, cfile, cal_dir, settings, default_color):
    """
    Calibrate the source tables of the images of a night.

    The sources in the images of each object are identified by their
    celestial coordinates (see :func:`_cross_match`), and calibrated with
    :func:`~pipnick.photometry.calibrate.iterate_colors`, so that the color
    of each source is computed from its calibrated magnitudes, averaged over
    the images in each band.  Sources without a measured color use the
    default color, and are flagged (``color_assumed``).

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
        The images; see :func:`_collect_images`.
    content : dict
        The calibration; see :func:`write_calibration`.
    cfile : :class:`pathlib.Path`
        The calibration file.
    cal_dir : :class:`pathlib.Path`
        The output directory.
    settings : dict
        The calibration settings; see :func:`calibrate_all`.
    default_color : dict
        The default colors, keyed by the color, e.g., ``{'B-V': 0.6}``.

    Returns
    -------
    list
        The paths to the calibrated tables.
    """
    calibrations = _calibrations(content)
    defaults = {parse_color(c): float(v) for c, v in default_color.items()}
    mag_type = content['calibration']['magnitude']
    chash = _file_hash(cfile)
    outputs = []
    usable = []
    for image in images:
        if image['band'] not in calibrations or image['airmass'] is None:
            log.warning(
                f'{image["path"].name}: its band is not calibrated; not calibrating its source '
                'table.'
            )
            continue
        usable.append(image)
    objects = sorted({str(image['phot'].meta.get('object')) for image in usable})
    for obj in objects:
        group = [image for image in usable if str(image['phot'].meta.get('object')) == obj]
        ids = _cross_match(group, settings['match_radius'])
        mags = []
        errs = []
        for image in group:
            m, e = _instrumental(image['phot'], mag_type)
            mags.append(m)
            errs.append(e)
        n = [len(m) for m in mags]
        mag, mag_err, color, assumed = iterate_colors(
            np.ma.concatenate(mags), np.ma.concatenate(errs),
            np.concatenate([np.full(k, float(i['airmass'])) for k, i in zip(n, group)]),
            np.concatenate(ids), np.concatenate([[i['band']] * k for k, i in zip(n, group)]),
            calibrations, default_color=defaults
        )
        start = 0
        for image, k, sid, m, e in zip(group, n, ids, mags, errs):
            phot = image['phot']
            cal = calibrations[image['band']]
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
                'image': image['path'].name, 'object': obj, 'band': image['band'],
                'color': cal['color_name'], 'magnitude': mag_type, 'system': 'Vega',
                'extinction': cal['extinction'], 'calibration': cfile.name,
                'calibration_hash': chash, 'airmass': image['airmass'],
                'errors': 'statistical only; fixed and fallback terms have no error'
            }
            ofile = cal_dir / f'{image["path"].stem}_cal.ecsv'
            tbl.write(ofile, format='ascii.ecsv', overwrite=True)
            outputs.append(ofile)
            log.info(
                f'Wrote {ofile.name}: {np.sum(np.logical_not(np.ma.getmaskarray(tbl["mag"])))} '
                f'calibrated sources, {np.sum(tbl["color_assumed"])} with an assumed color.'
            )
            start += k
    return outputs


def _summary_rows(night, content, status):
    """
    Return the summary rows of the bands of a calibration; see
    :func:`_summarize`.
    """
    rows = []
    for band, entry in content.get('band', {}).items():
        ok = entry.get('status') == 'ok'
        k = -1.
        if ok and 'airmass' in entry['terms']:
            k = entry['coefficients'][entry['terms'].index('airmass')]
        drift = entry.get('drift', [0., -1.])
        rows.append(
            {
                'night': night, 'band': band, 'status': entry.get('status', 'failed'),
                'n_stars': int(entry.get('n_stars', -1)), 'n_meas': int(entry.get('n_meas', -1)),
                'n_rejected': int(entry.get('n_rejected', -1)), 'rms': entry.get('rms', -1.),
                'rms_loo': entry.get('rms_loo', -1.), 'extinction': _extinction_source(entry),
                'k_airmass': k, 'drift': drift[0], 'drift_err': drift[1],
                'file': status, 'reason': entry.get('reason', entry.get('fallback_reason', '-'))
            }
        )
    return rows


def _summarize(summary, cal_dir):
    """
    Log the summary of the calibration of each band, write it to
    ``calibrate/summary.ecsv``, and end with a banner listing the bands whose
    extinction was not measured.
    """
    if len(summary) == 0:
        log.warning('No bands to calibrate.')
        return
    tbl = Table(rows=summary)
    for name in ['n_stars', 'n_meas', 'n_rejected', 'rms', 'rms_loo', 'k_airmass']:
        tbl[name] = np.ma.MaskedArray(tbl[name], mask=np.asarray(tbl[name]) < 0)
    tbl['drift'] = np.ma.MaskedArray(tbl['drift'], mask=np.asarray(tbl['drift_err']) < 0)
    tbl['drift_err'] = np.ma.MaskedArray(tbl['drift_err'], mask=np.asarray(tbl['drift_err']) < 0)
    for name in ['rms', 'rms_loo', 'k_airmass', 'drift', 'drift_err']:
        tbl[name].info.format = '.3f'
    cal_dir.mkdir(parents=True, exist_ok=True)
    tbl.write(cal_dir / 'summary.ecsv', format='ascii.ecsv', overwrite=True)
    columns = [
        'night', 'band', 'status', 'n_stars', 'n_meas', 'n_rejected', 'rms', 'rms_loo',
        'extinction', 'k_airmass', 'drift', 'file'
    ]
    log.info(
        'Calibration summary (also in calibrate/summary.ecsv):\n'
        + '\n'.join(tbl[columns].pformat(max_lines=-1, max_width=-1))
    )
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


def _write_qa(meas, content, ofile):
    """
    Write the QA figure of the calibration of a night: the residuals of each
    calibrated band against airmass, color, catalog magnitude, and time.

    Used measurements are filled circles, rejected measurements are crosses,
    and unusable measurements are not shown.  Bands whose extinction was not
    measured are labeled.
    """
    bands = [b for b, e in content['band'].items() if e.get('status') == 'ok']
    if len(bands) == 0:
        return
    fig = Figure(figsize=(14, 3 * len(bands)))
    axes = fig.subplots(len(bands), 4, squeeze=False)
    mjd0 = np.min(meas['mjd']) if len(meas) > 0 else 0.
    for row, band in zip(axes, bands):
        entry = content['band'][band]
        sel = meas['band'] == band
        resid = np.ma.filled(meas['resid'], 0.)
        err = np.ma.filled(meas['resid_err'], 0.)
        xs = {
            'airmass': np.asarray(meas['airmass']),
            entry['color']: np.ma.filled(meas['color'], 0.),
            f'{band} (catalog)': np.ma.filled(meas['mag_cat'], 0.),
            'hours': (np.asarray(meas['mjd']) - mjd0) * 24.
        }
        used = sel & np.asarray(meas['used'])
        rejected = sel & np.asarray(meas['rejected'])
        for ax, (label, x) in zip(row, xs.items()):
            ax.errorbar(x[used], resid[used], yerr=err[used], fmt='o', color='C0', ms=4)
            ax.scatter(x[rejected], resid[rejected], marker='x', color='C3', s=30)
            ax.axhline(0., color='k', lw=0.5)
            ax.set_xlabel(label)
            ax.set_ylabel(f'{band} residual (mag)')
        title = f'{band}: rms {entry["rms"]:.3f} mag, {entry["n_stars"]} stars'
        if 'rms_loo' in entry:
            title += f', leave-one-out rms {entry["rms_loo"]:.3f}'
        row[0].set_title(title, loc='left', fontsize=9)
        if 'fallback' in entry['source']:
            for ax in row:
                ax.text(
                    0.5, 0.92, 'EXTINCTION NOT MEASURED: nominal value used',
                    transform=ax.transAxes, ha='center', va='top', color='C3', fontsize=8,
                    fontweight='bold'
                )
        if 'drift' in entry:
            row[3].set_title(
                f'drift {entry["drift"][0]:+.4f} +/- {entry["drift"][1]:.4f} mag/h',
                loc='right', fontsize=9
            )
    fig.suptitle(
        f'{content["calibration"]["setup"]} {content["calibration"]["night"]}: calibrated minus '
        'catalog magnitudes'
    )
    fig.tight_layout()
    ofile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ofile)
    log.info(f'Wrote {ofile}')
