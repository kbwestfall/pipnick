"""
Camera definitions for the Nickel Telescope.

Each camera is defined by a TOML configuration file in ``pipnick/data/cameras``.
The files are found automatically, so supporting a new camera only requires
adding its configuration file; see :func:`available_cameras`.  The camera used
to obtain a raw frame is identified from its header with
:func:`identify_camera`, or a camera can be selected by name with
:func:`get_camera`.

Configuration schema
--------------------

All geometry and bad-pixel coordinates are in *unbinned*, 1-based (FITS)
detector pixels, and pixel ranges are inclusive ``[first, last]`` pairs.  The
binned and windowed equivalents for a given frame are derived from its header.

``[camera]``
    ``name`` (used to select the camera), ``description``, ``pixel_size``
    (micron), ``plate_scale`` (arcsec per unbinned pixel), ``saturation``
    (ADU), and ``fwhm_guess`` (arcsec).

``[identify]``
    Header keywords mapped to regular expressions.  A header is from this
    camera if all of the keywords are present and match.

``[keywords]``
    Generic quantities mapped to *ordered lists* of header keywords; the first
    keyword present in a header is used.  Quantities are ``object``,
    ``frametype``, ``obstype``, ``exptime``, ``filter``, ``airmass``, ``ra``,
    ``dec``, ``date``, ``ccdsum``, ``xbin``, ``ybin``, ``readspeed``, ``gain``,
    ``readnoise``, ``datasec``, ``biassec``, ``window``, and ``temperature``.
    Listing a standard keyword first means it is used as soon as it appears
    in the headers, with the camera-specific keywords as the fallback.

``[geometry]``
    ``shape`` (the full unbinned readout as ``[rows, columns]``),
    ``data_columns``, ``data_rows``, and ``overscan_columns`` (the columns used
    to measure the overscan level).  Binned pixels that straddle a boundary
    are excluded from both sides.  Used only when the header does not provide
    ``datasec`` or ``biassec``.

``[[readout]]``
    One entry per readout setting.  ``value`` is matched against the header
    ``readspeed`` keyword; ``mode`` is the generic name of the setting (e.g.,
    ``'slow'`` or ``'fast'``); ``gain`` (e-/ADU) and ``readnoise`` (e-) are
    used unless the header provides them; ``placeholder`` flags values that
    have not been measured; and ``saturation`` optionally overrides the camera
    value.

``[bpm]``
    The archival bad-pixel mask: ``columns`` and ``rows`` (lists of inclusive
    ranges), ``regions`` (polygons given as lists of ``[x, y]`` vertices, with
    pixel centers at integer coordinates), and an optional ``file`` (an
    unbinned mask image in ``pipnick/data/cameras``).

``[quality]``
    Thresholds used to reject frames; see :func:`NickelCamera.assess`.
    ``max_temperature`` (C), ``bias_level`` (the ``[min, max]`` expected level
    of closed-shutter frames, in ADU), ``bias_exptime_max`` (s),
    ``min_illumination`` (the minimum level of a flat above the bias level, in
    ADU), and ``max_saturated_fraction``.

``[frametypes]``
    ``default`` frame type, a ``map`` from explicit header frame-type values
    (see the ``frametype`` keyword) to pipnick frame types, and an ordered list
    of ``rule`` tables.  The first matching rule sets the frame type.  A rule
    may test ``obstype`` and ``object`` (case-insensitive regular
    expressions) and ``exptime_min`` and ``exptime_max`` (inclusive bounds, in
    seconds).
"""

from dataclasses import dataclass
from functools import cache
import importlib.resources
from pathlib import Path
import re
import tomllib

import numpy as np
from astropy.io import fits

from pipnick.pkg.exceptions import PipnickDataModelError, PipnickParameterError
from pipnick.utils import masking

__all__ = [
    'NickelCamera', 'Readout', 'available_cameras', 'get_camera', 'identify_camera', 'get_header'
]


@dataclass(frozen=True)
class Readout:
    """
    Properties of the readout setting used for a frame.
    """

    mode: str
    """
    Generic name of the readout setting, or ``'unknown'`` if the setting is
    not defined by the camera configuration.
    """

    gain: float | None
    """
    Gain in e-/ADU.
    """

    readnoise: float | None
    """
    Read noise in e-.
    """

    saturation: float
    """
    Saturation level in ADU.
    """

    placeholder: bool
    """
    The gain and read noise are placeholders that have not been measured.
    """


class NickelCamera:
    """
    A Nickel camera, as defined by its configuration file.

    Parameters
    ----------
    config : dict
        The camera configuration; see the module documentation for the schema.
    source : str, :class:`pathlib.Path`, optional
        The file the configuration was read from, used in error messages and
        to find the optional bad-pixel mask image.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if the configuration is missing required content.
    """

    required = {
        'camera': ['name', 'plate_scale', 'saturation'],
        'identify': [],
        'keywords': [],
        'geometry': ['shape', 'data_columns', 'data_rows', 'overscan_columns'],
        'frametypes': ['default'],
    }
    """
    Required configuration sections and their required entries.
    """

    def __init__(self, config, source=None):
        self.config = config
        self.source = None if source is None else Path(source)
        self._validate()

    @classmethod
    def from_file(cls, path):
        """
        Instantiate a camera from a TOML configuration file.

        Parameters
        ----------
        path : str, :class:`pathlib.Path`
            The configuration file.

        Returns
        -------
        :class:`NickelCamera`
            The camera.
        """
        _path = Path(path)
        with open(_path, 'rb') as f:
            return cls(tomllib.load(f), source=_path)

    def _validate(self):
        """
        Check that the configuration has the required content.
        """
        src = '' if self.source is None else f' in {self.source}'
        for section, entries in self.required.items():
            if section not in self.config:
                raise PipnickParameterError(f'Camera configuration{src} has no [{section}].')
            missing = [e for e in entries if e not in self.config[section]]
            if len(missing) > 0:
                raise PipnickParameterError(
                    f'Camera configuration{src} is missing {missing} '
                    f'from [{section}].'
                )
        if len(self.config['identify']) == 0:
            raise PipnickParameterError(
                f'Camera configuration{src} must define at least one '
                'keyword in [identify].'
            )

    def __repr__(self):
        return f'{self.__class__.__name__}({self.name})'

    @property
    def name(self):
        """
        The name of the camera.
        """
        return self.config['camera']['name']

    @property
    def description(self):
        """
        A short description of the camera.
        """
        return self.config['camera'].get('description', self.name)

    @property
    def pixel_size(self):
        """
        The size of the detector pixels in microns.
        """
        return self.config['camera'].get('pixel_size')

    @property
    def plate_scale(self):
        """
        The plate scale in arcsec per unbinned pixel.
        """
        return self.config['camera']['plate_scale']

    @property
    def saturation(self):
        """
        The saturation level in ADU.
        """
        return self.config['camera']['saturation']

    @property
    def fwhm_guess(self):
        """
        A guess for the seeing FWHM in arcsec.
        """
        return self.config['camera'].get('fwhm_guess')

    @property
    def shape(self):
        """
        The shape of the full unbinned readout, as ``(rows, columns)``.
        """
        return tuple(self.config['geometry']['shape'])

    def matches(self, hdr):
        """
        Check if a header was produced by this camera.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The header to check.

        Returns
        -------
        bool
            True if all of the ``[identify]`` keywords are present in the
            header and match.
        """
        return all(
            key in hdr and re.search(pattern, str(hdr[key])) is not None for key,
            pattern in self.config['identify'].items()
        )

    def header_keyword(self, hdr, quantity):
        """
        Find the header keyword that provides a quantity.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The header to search.
        quantity : str
            The generic name of the quantity; see the ``[keywords]`` section of
            the configuration.

        Returns
        -------
        str, None
            The first keyword for the quantity that is present in the header,
            or None if there is none.
        """
        for key in self.config['keywords'].get(quantity, []):
            if key in hdr:
                return key
        return None

    def header_value(self, hdr, quantity, default=None):
        """
        Get the value of a quantity from a header.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The header to search.
        quantity : str
            The generic name of the quantity; see the ``[keywords]`` section of
            the configuration.
        default : object, optional
            The value to return if the header does not provide the quantity.

        Returns
        -------
        object
            The header value.
        """
        key = self.header_keyword(hdr, quantity)
        return default if key is None else hdr[key]

    def binning(self, hdr):
        """
        Get the binning of a frame.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The frame header.

        Returns
        -------
        tuple
            The binning along the columns (x) and rows (y), as positive
            integers.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
            Raised if the header does not provide the binning.
        """
        ccdsum = self.header_value(hdr, 'ccdsum')
        if ccdsum is not None:
            return tuple(int(b) for b in str(ccdsum).split())
        xbin = self.header_value(hdr, 'xbin')
        ybin = self.header_value(hdr, 'ybin')
        if xbin is None or ybin is None:
            raise PipnickDataModelError(f'Unable to determine the binning for {self.name} data.')
        return abs(int(xbin)), abs(int(ybin))

    def window(self, hdr):
        """
        Get the readout window of a frame.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The frame header.

        Returns
        -------
        tuple
            The inclusive, 1-based, unbinned first and last columns and rows of
            the readout window, as ``(x1, x2, y1, y2)``.  If the header does not
            define a window, the full readout is returned.
        """
        window = self.header_value(hdr, 'window')
        if window is None:
            ny, nx = self.shape
            return 1, nx, 1, ny
        return tuple(int(w) for w in str(window).split())

    def is_windowed(self, hdr):
        """
        Check if a frame is a subframe of the full readout.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The frame header.

        Returns
        -------
        bool
            True if the readout window is smaller than the full readout.
        """
        ny, nx = self.shape
        return self.window(hdr) != (1, nx, 1, ny)

    def readout(self, hdr):
        """
        Get the readout properties of a frame.

        Values provided by the header (see the ``gain`` and ``readnoise``
        keywords) take precedence over the configuration.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The frame header.

        Returns
        -------
        :class:`Readout`
            The readout properties.  If the readout setting is not defined in
            the configuration, the mode is ``'unknown'``.
        """
        value = self.header_value(hdr, 'readspeed')
        entry = {}
        if value is not None:
            for _entry in self.config.get('readout', []):
                if str(value).strip() == str(_entry['value']):
                    entry = _entry
                    break
        gain = self.header_value(hdr, 'gain', default=entry.get('gain'))
        readnoise = self.header_value(hdr, 'readnoise', default=entry.get('readnoise'))
        # The values are only placeholders if they come from the configuration
        placeholder = (
            entry.get('placeholder', False)
            and self.header_keyword(hdr, 'gain') is None
            and self.header_keyword(hdr, 'readnoise') is None
        )
        return Readout(
            mode=entry.get('mode', 'unknown'), gain=gain, readnoise=readnoise,
            saturation=entry.get('saturation', self.saturation), placeholder=placeholder
        )

    @staticmethod
    def _binned_range(rng, start, nbin, npix):
        """
        Convert an unbinned pixel range to the binned pixels fully within it.

        Parameters
        ----------
        rng : list
            Inclusive, 1-based ``[first, last]`` unbinned pixel range.
        start : int
            The first (1-based) unbinned pixel of the readout window.
        nbin : int
            The binning.
        npix : int
            The number of binned pixels in the frame.

        Returns
        -------
        tuple, None
            The inclusive, 1-based ``(first, last)`` binned pixel range, or None
            if no binned pixel falls fully within the range.
        """
        # Binned pixel j (1-based) covers the unbinned pixels
        # start + (j-1)*nbin to start + j*nbin - 1
        first = max(int(np.ceil((rng[0] - start) / nbin + 1)), 1)
        last = min(int(np.floor((rng[1] - start + 1) / nbin)), npix)
        return None if first > last else (first, last)

    def raw_shape(self, hdr):
        """
        Get the expected shape of a raw frame.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The frame header.

        Returns
        -------
        tuple
            The binned shape of the readout window, as ``(rows, columns)``.
        """
        xbin, ybin = self.binning(hdr)
        x1, x2, y1, y2 = self.window(hdr)
        return (y2 - y1 + 1) // ybin, (x2 - x1 + 1) // xbin

    def _check_shape(self, hdr):
        """
        Check that the header image size matches the expected frame shape.
        """
        shape = self.raw_shape(hdr)
        if 'NAXIS1' in hdr and 'NAXIS2' in hdr and (hdr['NAXIS2'], hdr['NAXIS1']) != shape:
            raise PipnickDataModelError(
                f'{self.name} frame shape ({hdr["NAXIS2"]}, {hdr["NAXIS1"]}) does not match '
                f'the shape expected from its binning and readout window: {shape}.'
            )
        return shape

    def sections(self, hdr):
        """
        Get the data and overscan sections of a raw frame.

        Sections provided by the header (see the ``datasec`` and ``biassec``
        keywords) take precedence over the configured geometry.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The frame header.

        Returns
        -------
        datasec : str
            The data section, as a FITS section string (e.g.,
            ``'[1:1024,1:1024]'``).
        biassec : str, None
            The overscan section used to measure the overscan level, or None if
            the frame has no overscan (e.g., a readout window that excludes
            it).

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
            Raised if the frame shape does not match its binning and readout
            window, or if the frame contains no data pixels.
        """
        datasec = self.header_value(hdr, 'datasec')
        biassec = self.header_value(hdr, 'biassec')
        if datasec is not None and biassec is not None:
            return datasec, biassec

        ny, nx = self._check_shape(hdr)
        xbin, ybin = self.binning(hdr)
        x1, x2, y1, y2 = self.window(hdr)
        geom = self.config['geometry']
        rows = self._binned_range(geom['data_rows'], y1, ybin, ny)
        if datasec is None:
            cols = self._binned_range(geom['data_columns'], x1, xbin, nx)
            if cols is None or rows is None:
                raise PipnickDataModelError(f'{self.name} frame contains no data pixels.')
            datasec = f'[{cols[0]}:{cols[1]},{rows[0]}:{rows[1]}]'
        if biassec is None:
            cols = self._binned_range(geom['overscan_columns'], x1, xbin, nx)
            if cols is not None and rows is not None:
                biassec = f'[{cols[0]}:{cols[1]},{rows[0]}:{rows[1]}]'
        return datasec, biassec

    @staticmethod
    def _range_indices(ranges):
        """
        Convert a list of inclusive, 1-based ranges to 0-based indices.
        """
        if len(ranges) == 0:
            return None
        return np.concatenate([np.arange(r[0] - 1, r[1]) for r in ranges])

    @cache
    def unbinned_bpm(self):
        """
        Construct the archival bad-pixel mask for the full unbinned readout.

        Returns
        -------
        :class:`numpy.ndarray`
            Boolean mask with the shape of the full unbinned readout (see
            :attr:`shape`); bad pixels are True.
        """
        bpm = self.config.get('bpm', {})
        # Convert from 1-based to 0-based pixel coordinates
        regions = [np.asarray(r, dtype=float) - 1 for r in bpm.get('regions', [])]
        mask = masking.create_mask(
            self.shape, bad_cols=self._range_indices(bpm.get('columns', [])),
            bad_regions=regions if len(regions) > 0 else None
        )
        rows = self._range_indices(bpm.get('rows', []))
        if rows is not None:
            mask[rows] = True
        if bpm.get('file', '') != '':
            mask_file = (
                (self.source.parent if self.source is not None else Path('.')) / bpm['file']
            )
            file_mask = fits.getdata(mask_file).astype(bool)
            if file_mask.shape != self.shape:
                raise PipnickParameterError(
                    f'Bad-pixel mask in {mask_file} has shape '
                    f'{file_mask.shape}; expected {self.shape}.'
                )
            mask |= file_mask
        return mask

    def bpm(self, hdr):
        """
        Construct the bad-pixel mask for a raw frame.

        The unbinned mask is cut to the frame's readout window and binned; a
        binned pixel is bad if any of the unbinned pixels it contains is bad.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The frame header.

        Returns
        -------
        :class:`numpy.ndarray`
            Boolean mask with the shape of the raw frame; bad pixels are True.
        """
        ny, nx = self._check_shape(hdr)
        xbin, ybin = self.binning(hdr)
        x1, x2, y1, y2 = self.window(hdr)
        mask = self.unbinned_bpm()[y1 - 1:y1 - 1 + ny * ybin, x1 - 1:x1 - 1 + nx * xbin]
        return mask.reshape(ny, ybin, nx, xbin).any(axis=(1, 3))

    def classify(self, hdr):
        """
        Determine the frame type from a header.

        An explicit frame type in the header (see the ``frametype`` keyword) is
        used if it is defined in ``[frametypes.map]``.  Otherwise, the
        ``[frametypes]`` rules are tested in order, and the first that matches
        sets the frame type.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The frame header.

        Returns
        -------
        str
            The frame type.
        """
        frametypes = self.config['frametypes']
        explicit = self.header_value(hdr, 'frametype')
        if explicit is not None:
            mapped = frametypes.get('map', {}).get(str(explicit).strip().lower())
            if mapped is not None:
                return mapped

        obstype = self.header_value(hdr, 'obstype')
        obj = self.header_value(hdr, 'object')
        exptime = self.header_value(hdr, 'exptime')
        for rule in frametypes.get('rule', []):
            if 'obstype' in rule and (
                obstype is None
                or re.search(rule['obstype'], str(obstype).strip(), flags=re.IGNORECASE) is None
            ):
                continue
            if 'object' in rule and (
                obj is None
                or re.search(rule['object'], str(obj).strip(), flags=re.IGNORECASE) is None
            ):
                continue
            if 'exptime_min' in rule and (exptime is None or exptime < rule['exptime_min']):
                continue
            if 'exptime_max' in rule and (exptime is None or exptime > rule['exptime_max']):
                continue
            return rule['frametype']
        return frametypes['default']

    @staticmethod
    def _section_slice(section):
        """
        Convert a FITS section string to a numpy slice.
        """
        (x1, x2), (y1, y2) = [map(int, s.split(':')) for s in section.strip('[]').split(',')]
        return np.s_[y1 - 1:y2, x1 - 1:x2]

    def assess(self, hdr, frametype, data=None):
        """
        Check the quality of a frame.

        The checks use the thresholds in the ``[quality]`` section of the
        configuration.  Using the header, frames are rejected if the detector
        was too warm or the readout setting is not defined.  If the image data
        are provided, frames are also rejected if they are empty (all zeros,
        e.g., from a failed readout), have too many saturated pixels, are
        closed-shutter frames (biases and darks) with a level above the bias
        level, or are flats with too little illumination.  In addition, frames
        that are labeled as calibrations but were exposed (frame type
        ``'unknown'``) are reclassified as biases if their exposure time and
        level are consistent with a bias.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The frame header.
        frametype : str
            The frame type; see :func:`classify`.
        data : :class:`numpy.ndarray`, optional
            The raw image data.  If None, only the header checks are performed.

        Returns
        -------
        frametype : str
            The frame type, which is changed only if the frame is reclassified.
        reject : bool
            True if the frame should be excluded from the reduction.
        notes : list
            Strings describing why the frame was rejected or reclassified.
        """
        quality = self.config.get('quality', {})
        notes = []
        reject = False

        temperature = self._to_float(self.header_value(hdr, 'temperature'))
        if (
            'max_temperature' in quality and temperature is not None
            and temperature > quality['max_temperature']
        ):
            reject = True
            notes.append(f'detector temperature {temperature:.1f} C is too warm')

        readout = self.readout(hdr)
        if readout.mode == 'unknown':
            reject = True
            notes.append(f'readout setting {self.header_value(hdr, "readspeed")} is undefined')

        if data is None:
            return frametype, reject, notes

        if np.all(data == 0):
            notes.append('empty image')
            return frametype, True, notes

        datasec, _ = self.sections(hdr)
        data = data[self._section_slice(datasec)]
        level = float(np.median(data))
        if 'max_saturated_fraction' in quality:
            frac = float(np.mean(data >= readout.saturation))
            if frac > quality['max_saturated_fraction']:
                reject = True
                notes.append(f'{100 * frac:.0f}% of pixels are saturated')

        if 'bias_level' not in quality:
            return frametype, reject, notes
        bias_lo, bias_hi = quality['bias_level']
        exptime = self._to_float(self.header_value(hdr, 'exptime'))
        if (
            frametype == 'unknown' and bias_lo <= level <= bias_hi and exptime is not None
            and exptime <= quality.get('bias_exptime_max', 0.)
        ):
            frametype = 'bias'
            notes.append(f'reclassified as bias (level {level:.0f} ADU)')
        if frametype in ['bias', 'dark'] and level > bias_hi:
            reject = True
            notes.append(f'level {level:.0f} ADU is too high for a {frametype}')
        if (
            frametype in ['domeflat', 'skyflat']
            and level < bias_hi + quality.get('min_illumination', 0.)
        ):
            reject = True
            notes.append(f'level {level:.0f} ADU is too low for a flat')
        return frametype, reject, notes

    def window_label(self, hdr):
        """
        Construct a label for the readout window of a frame.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The frame header.

        Returns
        -------
        str
            ``'full'`` for the full readout, or the unbinned window as a FITS
            section string (e.g., ``'[960:1205,840:1069]'``).
        """
        if not self.is_windowed(hdr):
            return 'full'
        x1, x2, y1, y2 = self.window(hdr)
        return f'[{x1}:{x2},{y1}:{y2}]'

    @staticmethod
    def _to_float(value):
        """
        Convert a header value to a float, returning None if not possible.
        """
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _to_str(value):
        """
        Convert a header value to a stripped string, returning None for None.
        """
        return None if value is None else str(value).strip()

    def metadata(self, hdr):
        """
        Collect the metadata for a frame, normalized across cameras.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The frame header.

        Returns
        -------
        dict
            The frame metadata: ``frametype``, ``object``, ``filter``,
            ``exptime`` (s), ``airmass``, ``ra``, ``dec``, ``date`` (ISO
            format), ``xbin``, ``ybin`` (positive integers), ``readmode``,
            ``naxis1``, ``naxis2``, and ``window`` (see
            :func:`window_label`).  Quantities missing from the header are
            None.
        """
        xbin, ybin = self.binning(hdr)
        return {
            'frametype': self.classify(hdr),
            'object': self._to_str(self.header_value(hdr, 'object')),
            'filter': self._to_str(self.header_value(hdr, 'filter')),
            'exptime': self._to_float(self.header_value(hdr, 'exptime')),
            'airmass': self._to_float(self.header_value(hdr, 'airmass')),
            'ra': self._to_str(self.header_value(hdr, 'ra')),
            'dec': self._to_str(self.header_value(hdr, 'dec')),
            'date': self._to_str(self.header_value(hdr, 'date')),
            'xbin': xbin,
            'ybin': ybin,
            'readmode': self.readout(hdr).mode,
            'naxis1': hdr.get('NAXIS1'),
            'naxis2': hdr.get('NAXIS2'),
            'window': self.window_label(hdr),
        }


def _config_files():
    """
    Find the camera configuration files distributed with pipnick.
    """
    root = importlib.resources.files('pipnick') / 'data' / 'cameras'
    return sorted(Path(str(f)) for f in root.iterdir() if f.name.endswith('.toml'))


@cache
def _registry():
    """
    Load all of the cameras distributed with pipnick, keyed by name.
    """
    registry = {}
    for f in _config_files():
        camera = NickelCamera.from_file(f)
        if camera.name in registry:
            raise PipnickParameterError(f'Camera name {camera.name} in {f} is not unique.')
        registry[camera.name] = camera
    return registry


def available_cameras():
    """
    List the cameras defined by pipnick.

    Returns
    -------
    list
        The names of the available cameras.
    """
    return list(_registry().keys())


def get_camera(name):
    """
    Get a camera by name.

    Parameters
    ----------
    name : str
        The name of the camera; see :func:`available_cameras`.

    Returns
    -------
    :class:`NickelCamera`
        The camera.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if the camera is not defined.
    """
    registry = _registry()
    if name not in registry:
        raise PipnickParameterError(
            f'Unknown camera {name}.  Available cameras are: '
            f'{", ".join(registry.keys())}.'
        )
    return registry[name]


def identify_camera(inp):
    """
    Identify the camera that produced a raw frame.

    Parameters
    ----------
    inp : str, :class:`pathlib.Path`, :class:`astropy.io.fits.Header`
        The raw file, or its header.

    Returns
    -------
    :class:`NickelCamera`
        The camera.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if no camera, or more than one camera, matches the header.
    """
    hdr = get_header(inp)
    matched = [camera for camera in _registry().values() if camera.matches(hdr)]
    if len(matched) == 1:
        return matched[0]
    if len(matched) == 0:
        raise PipnickDataModelError('Unable to determine the camera used to obtain the data.')
    raise PipnickDataModelError(
        'Header matches more than one camera: '
        f'{", ".join(c.name for c in matched)}.'
    )


def get_header(inp):
    """
    Convenience function used to get a header from the input.

    Parameters
    ----------
    inp : str, :class:`pathlib.Path`, :class:`astropy.io.fits.Header`
        File or header to read.  If a Header object, the object is simply
        returned.

    Returns
    -------
    hdr : :class:`astropy.io.fits.Header`
        Fits header
    """
    if isinstance(inp, fits.Header):
        return inp
    return fits.getheader(Path(inp).absolute())
