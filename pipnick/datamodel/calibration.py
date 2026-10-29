"""
The photometric calibration of a night.

The calibration of each night is written by
:func:`~pipnick.pipelines.calibrate.calibrate_all` to a TOML file,
``calibrate/<setup>_<night>.toml``, with a ``[calibration]`` table (the setup,
night, catalog, magnitude, creation time, pipnick version, and the hash of the
inputs) and a ``[band.<band>]`` table for each band.  :class:`Calibration`
reads and writes it.  The coefficients (and the terms) may be edited by the
user; the ``fit_hash`` entry, the hash of the content as written, is used to
detect edits (see :func:`calibration_hash`).
"""

import json
import tomllib
from pathlib import Path

import numpy as np

from pipnick.photometry.calibrate import check_terms, zero_point
from pipnick.photometry.landolt import parse_color
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickParameterError
from pipnick.pkg.provenance import file_hash, json_hash


def toml_value(value):
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
        return '[' + ', '.join(toml_value(v) for v in value) + ']'
    raise PipnickParameterError(f'Cannot write {type(value)} to the calibration file.')


def _toml_text(content):
    """
    Construct the text of a calibration file; see :meth:`Calibration.write`.
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
    lines += [f'{k} = {toml_value(v)}' for k, v in cal.items() if k != 'fit_hash']
    if 'fit_hash' in cal:
        lines.append(f'fit_hash = {toml_value(cal["fit_hash"])}')
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
        lines += [f'{k} = {toml_value(v)}' for k, v in entry.items()]
    return '\n'.join(lines) + '\n'


def calibration_hash(content):
    """
    Return the hash of the content of a calibration file, excluding its
    ``fit_hash`` entry.
    """
    _content = json.loads(json.dumps(content))
    _content.get('calibration', {}).pop('fit_hash', None)
    return json_hash(_content)


def header_calibration_state(header, cfile):
    """
    Determine if the photometric calibration recorded in the header of an
    image is up to date with the calibration file.

    This is the check of the application of the calibration (see
    :func:`~pipnick.photometry.apply.apply_to_headers`) and of the status of
    the data set (see :mod:`pipnick.report.status`).

    Parameters
    ----------
    header : :class:`astropy.io.fits.Header`
        The header of the image.
    cfile : :class:`pathlib.Path`
        The calibration file of its night.

    Returns
    -------
    str
        ``'ok'``, ``'stale'`` (the calibration file has changed since it was
        recorded, ``PNCALHSH``), or ``'missing'`` (no calibration is
        recorded).
    """
    if 'PNCALHSH' not in header:
        return 'missing'
    return 'ok' if header['PNCALHSH'] == file_hash(cfile) else 'stale'


def extinction_source(entry):
    """
    Return how the extinction of a band was set: ``'fit'``, ``'fixed'``,
    ``'fallback'``, ``'none'`` (no airmass term), or ``'failed'``.
    """
    if entry.get('status') != 'ok':
        return 'failed'
    if 'airmass' not in entry['terms']:
        return 'none'
    return entry['source'][entry['terms'].index('airmass')]


class Calibration:
    """
    The photometric calibration of a night.

    Parameters
    ----------
    content : dict
        The calibration: the ``calibration`` and ``band`` tables.
    path : :class:`pathlib.Path`, optional
        The calibration file.
    edited : bool, optional
        True if the content of the file no longer matches its ``fit_hash``,
        i.e., it has been edited.
    """
    def __init__(self, content, path=None, edited=False):
        self.content = content
        self.path = path
        self.edited = edited

    @classmethod
    def read(cls, path):
        """
        Read a calibration file.

        Parameters
        ----------
        path : :class:`pathlib.Path`
            The file.

        Returns
        -------
        :class:`Calibration`
            The calibration; see :attr:`Calibration.edited`.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
            Raised if the file cannot be parsed, or the coefficients of a band
            do not match its terms.
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
                raise PipnickDataModelError(
                    f'{path}: the covariance of band {band} is not {n}x{n}.'
                )
        edited = calibration_hash(content) != content['calibration'].get('fit_hash')
        return cls(content, path=Path(path), edited=edited)

    def write(self, path=None):
        """
        Write the calibration file, with its ``fit_hash`` (see
        :func:`calibration_hash`).

        Parameters
        ----------
        path : :class:`pathlib.Path`, optional
            The file; if None, the file it was read from.
        """
        _path = self.path if path is None else Path(path)
        self.content['calibration'].pop('fit_hash', None)
        # The hash is computed from the content as it is read back, so that it
        # does not depend on how the values are formatted
        parsed = tomllib.loads(_toml_text(self.content))
        self.content['calibration']['fit_hash'] = calibration_hash(parsed)
        _path.parent.mkdir(parents=True, exist_ok=True)
        _path.write_text(_toml_text(self.content))
        self.path = _path
        self.edited = False

    def bands(self):
        """
        Return the calibrations of the calibrated bands, in the form used by
        :func:`~pipnick.photometry.calibrate.iterate_colors`.

        Returns
        -------
        dict
            For each band: its terms (``terms``), coefficients (``coeff``),
            their covariance (``cov``), its color (``color``, parsed, and
            ``color_name``), and how its extinction was set
            (``extinction``).
        """
        calibrations = {}
        for band, entry in self.content.get('band', {}).items():
            if entry.get('status') != 'ok':
                continue
            n = len(entry['terms'])
            calibrations[band] = {
                'terms': entry['terms'], 'coeff': np.asarray(entry['coefficients'], dtype=float),
                'cov': np.asarray(entry.get('covariance', np.zeros((n, n))), dtype=float),
                'color': parse_color(entry['color']), 'color_name': entry['color'],
                'extinction': extinction_source(entry)
            }
        return calibrations

    def summary_rows(self, night, action):
        """
        Return the summary rows of the bands; see
        :func:`~pipnick.pipelines.calibrate.calibrate_all`.

        Parameters
        ----------
        night : str
            The night.
        action : str
            What was done with the calibration file, e.g., ``'fit'``.

        Returns
        -------
        list
            One dictionary for each band.  The zero point (``zp``, with its
            error ``zp_err``) is that above the atmosphere (airmass 0) for
            color 0; see :func:`~pipnick.photometry.calibrate.zero_point`.
        """
        rows = []
        for band, entry in self.content.get('band', {}).items():
            ok = entry.get('status') == 'ok'
            k = None
            zp = {'zp': None, 'zp_err': None}
            if ok:
                terms = entry['terms']
                if 'airmass' in terms:
                    k = entry['coefficients'][terms.index('airmass')]
                n = len(terms)
                zp = zero_point(
                    terms, entry['coefficients'], entry.get('covariance', np.zeros((n, n))), 0.
                )
            drift = entry.get('drift', [None, None])
            rows.append(
                {
                    'night': night, 'band': band, 'status': entry.get('status', 'failed'),
                    'n_stars': entry.get('n_stars'), 'n_meas': entry.get('n_meas'),
                    'n_rejected': entry.get('n_rejected'), 'rms': entry.get('rms'),
                    'rms_loo': entry.get('rms_loo'), 'extinction': extinction_source(entry),
                    'k_airmass': k, 'zp': zp['zp'], 'zp_err': zp['zp_err'], 'drift': drift[0],
                    'drift_err': drift[1], 'action': action,
                    'reason': entry.get('reason', entry.get('fallback_reason', '-'))
                }
            )
        return rows
