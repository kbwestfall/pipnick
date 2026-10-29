"""
The settings of the processing steps.

Each step that has settings collects them in a frozen dataclass:
:class:`AstrometrySettings`, :class:`PhotometrySettings`,
:class:`CalibrationSettings`, and :class:`MosaicSettings`.  The values are
normalized (e.g., converted to ``float``) and checked (see ``validate``) when
the settings are constructed, so invalid settings are rejected before
anything is processed.  A change to a setting is made with
:func:`dataclasses.replace`, which constructs, and so checks, new settings.

The ``hash`` of the settings of the photometry, the calibration, and the
mosaics is recorded with their products, and decides whether a product is up
to date.  It is the :func:`~pipnick.pkg.provenance.json_hash` of the fields
(see ``as_dict``), so the names, values, and types of
the fields are part of the on-disk record of each product, and must not
change.
"""

import dataclasses

from pipnick.photometry.calibrate import check_terms
from pipnick.photometry.landolt import default_colors, parse_color
from pipnick.pkg.exceptions import PipnickParameterError
from pipnick.pkg.provenance import json_hash


def _check_option(name, value, options):
    """
    Check that a setting is one of its options.
    """
    if value not in options:
        raise PipnickParameterError(f'Invalid {name} {value}; options are {options}.')


def _check_positive(name, value):
    """
    Check that a setting is positive.
    """
    if value <= 0:
        raise PipnickParameterError(f'Invalid {name} {value}; must be positive.')


def _resolve_fwhm(settings, camera):
    """
    Return settings with the FWHM set to the typical value for the camera, if
    it is not set.
    """
    if settings.fwhm is not None:
        return settings
    return dataclasses.replace(settings, fwhm=camera.fwhm_guess)


class _Settings:
    """
    The methods shared by the settings dataclasses.
    """
    def as_dict(self):
        """
        Return the settings as a dictionary, e.g., as recorded in the metadata
        of the source tables.
        """
        return dataclasses.asdict(self)

    def hash(self):
        """
        Return the hash of the settings; see
        :func:`~pipnick.pkg.provenance.json_hash`.
        """
        return json_hash(self.as_dict())


@dataclasses.dataclass(frozen=True)
class AstrometrySettings(_Settings):
    """
    The settings of the astrometric calibration; see
    :func:`~pipnick.pipelines.astrometry.astrometry_all` for their meaning.

    The settings are not recorded with the solutions; the detected sources
    and the astrometry.net hints set the key of the cached astrometry.net
    results (see :func:`~pipnick.astrometry.cache.anet_input_hash`).
    """
    max_sources: int = 100
    detect_thresh: float = 5.
    fwhm: float = None
    upload_image: bool = False
    max_submissions: int = None
    group_radius: float = 1.
    min_matches: int = 4
    max_rms: float = 0.3
    search: float = 60.
    wide_search: float = 10.
    wide_min_matches: int = 10

    def __post_init__(self):
        for name in ['max_sources', 'min_matches', 'wide_min_matches']:
            object.__setattr__(self, name, int(getattr(self, name)))
        for name in ['detect_thresh', 'group_radius', 'max_rms', 'search', 'wide_search']:
            object.__setattr__(self, name, float(getattr(self, name)))
        if self.fwhm is not None:
            object.__setattr__(self, 'fwhm', float(self.fwhm))
        if self.max_submissions is not None:
            object.__setattr__(self, 'max_submissions', int(self.max_submissions))
        object.__setattr__(self, 'upload_image', bool(self.upload_image))
        self.validate()

    def validate(self):
        """
        Check the settings.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
            Raised if a setting is invalid.
        """
        for name in ['max_sources', 'detect_thresh', 'group_radius', 'max_rms', 'search']:
            _check_positive(name, getattr(self, name))
        if self.fwhm is not None:
            _check_positive('fwhm', self.fwhm)
        if self.wide_search < 0:
            raise PipnickParameterError(
                f'Invalid wide_search {self.wide_search}; must not be negative.'
            )
        if self.max_submissions is not None and self.max_submissions < 0:
            raise PipnickParameterError(
                f'Invalid max_submissions {self.max_submissions}; must not be negative.'
            )

    def resolve(self, camera):
        """
        Return the settings with the FWHM set to the typical value for the
        camera (``fwhm_guess``), if it is not set.

        Parameters
        ----------
        camera : :class:`~pipnick.cameras.cameras.NickelCamera`
            The camera.

        Returns
        -------
        :class:`AstrometrySettings`
            The settings.
        """
        return _resolve_fwhm(self, camera)


@dataclasses.dataclass(frozen=True)
class PhotometrySettings(_Settings):
    """
    The settings of the photometry; see
    :func:`~pipnick.pipelines.photometry.photometry_all` for their meaning.

    The FWHM defaults to the typical value for the camera; see
    :meth:`resolve`.  The settings, with the FWHM resolved, and their hash
    are recorded in the metadata of each source table.
    """
    thresh: float = 8.
    fwhm: float = None
    fittype: str = 'circ'
    mode: str = 'all'
    stack: str = 'mean'
    group: bool = False
    min_sep: float = 1.
    group_sep: float = 2.
    aper_radius: float = 2.
    psf_thresh: float = 15.

    def __post_init__(self):
        for name in ['thresh', 'min_sep', 'group_sep', 'aper_radius', 'psf_thresh']:
            object.__setattr__(self, name, float(getattr(self, name)))
        if self.fwhm is not None:
            object.__setattr__(self, 'fwhm', float(self.fwhm))
        object.__setattr__(self, 'group', bool(self.group))
        self.validate()

    def validate(self):
        """
        Check the settings.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
            Raised if a setting is invalid.
        """
        _check_option('fittype', self.fittype, ['circ', 'ellip'])
        _check_option('mode', self.mode, ['all', 'new'])
        _check_option('stack', self.stack, ['mean', 'median'])
        for name in ['thresh', 'min_sep', 'group_sep', 'aper_radius', 'psf_thresh']:
            _check_positive(name, getattr(self, name))
        if self.fwhm is not None:
            _check_positive('fwhm', self.fwhm)

    def resolve(self, camera):
        """
        Return the settings with the FWHM set to the typical value for the
        camera (``fwhm_guess``), if it is not set.

        Parameters
        ----------
        camera : :class:`~pipnick.cameras.cameras.NickelCamera`
            The camera.

        Returns
        -------
        :class:`PhotometrySettings`
            The settings.
        """
        return _resolve_fwhm(self, camera)

    def hash(self):
        """
        Return the hash of the settings, recorded in the metadata of the
        source tables (``settings_hash``); see
        :func:`~pipnick.pkg.provenance.json_hash`.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
            Raised if the FWHM has not been resolved; see :meth:`resolve`.
        """
        if self.fwhm is None:
            raise PipnickParameterError('Resolve the FWHM before hashing the settings.')
        return super().hash()


@dataclasses.dataclass(frozen=True)
class CalibrationSettings(_Settings):
    """
    The settings of the photometric calibration; see
    :func:`~pipnick.pipelines.calibrate.calibrate_all` for their meaning.

    The terms are put in their standard order (see
    :func:`~pipnick.photometry.calibrate.check_terms`), and the colors given
    are combined with the default colors (see
    :func:`~pipnick.photometry.landolt.default_colors`).  The fixed terms and
    the default colors of the sources are not settings of the fit; they are
    hashed with the settings into the input hash of each night (see
    :func:`~pipnick.pipelines.calibrate.calibrate_all`).
    """
    terms: list = dataclasses.field(default_factory=lambda: ['zp', 'airmass', 'color'])
    colors: dict = None
    mag: str = 'psf'
    match_radius: float = 2.
    max_neighbor_flux: float = 0.01
    residual_radius: float = 1.5
    max_cat_err: float = 0.05
    min_cat_obs: int = 1
    sys_err: float = 0.01
    rej: float = 3.
    maxiter: int = 5
    min_stars: int = 3
    min_airmass_span: float = 0.2
    min_color_span: float = 0.3
    force: bool = False
    fallback: bool = True
    drift_warn: float = 0.02

    def __post_init__(self):
        object.__setattr__(self, 'terms', check_terms(list(self.terms)))
        colors = default_colors()
        colors.update({} if self.colors is None else self.colors)
        object.__setattr__(self, 'colors', colors)
        for name in [
            'match_radius', 'max_neighbor_flux', 'residual_radius', 'max_cat_err', 'sys_err',
            'rej', 'min_airmass_span', 'min_color_span', 'drift_warn'
        ]:
            object.__setattr__(self, name, float(getattr(self, name)))
        for name in ['min_cat_obs', 'maxiter', 'min_stars']:
            object.__setattr__(self, name, int(getattr(self, name)))
        for name in ['force', 'fallback']:
            object.__setattr__(self, name, bool(getattr(self, name)))
        self.validate()

    def validate(self):
        """
        Check the settings.

        The colors are parsed (see
        :func:`~pipnick.photometry.landolt.parse_color`), but a color that does
        not include its band is accepted; see
        :func:`~pipnick.pipelines.calibrate.calibrate_all`.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
            Raised if a setting is invalid.
        """
        _check_option('mag', self.mag, ['psf', 'aper'])
        for color in self.colors.values():
            parse_color(color)
        for name in ['match_radius', 'residual_radius', 'rej']:
            _check_positive(name, getattr(self, name))
        if self.maxiter < 0:
            raise PipnickParameterError(f'Invalid maxiter {self.maxiter}; must not be negative.')


@dataclasses.dataclass(frozen=True)
class MosaicSettings(_Settings):
    """
    The settings of the mosaics; see
    :func:`~pipnick.pipelines.mosaic.mosaic_all` for their meaning.

    The settings are hashed into the hash of each mosaic (``PNMOSHSH``).  The
    reference image and the pixel scale are hashed as they are given.
    """
    psf_match: str = 'mosaic'
    kernel: str = 'moffat'
    psf_margin: float = 1.05
    scale: str = 'sources'
    reference: str = None
    match_radius: float = 1.
    min_common: int = 3
    min_snr: float = 20.
    background: str = 'constant'
    match_background: bool = False
    reproject: str = 'adaptive'
    pixel_scale: float = None
    rej: float = 3.
    rej_frac: float = 0.1

    def __post_init__(self):
        for name in ['psf_margin', 'match_radius', 'min_snr', 'rej_frac']:
            object.__setattr__(self, name, float(getattr(self, name)))
        object.__setattr__(self, 'min_common', int(self.min_common))
        object.__setattr__(self, 'match_background', bool(self.match_background))
        if self.rej is not None:
            object.__setattr__(self, 'rej', float(self.rej))
        self.validate()

    def validate(self):
        """
        Check the settings.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
            Raised if a setting is invalid.
        """
        _check_option('psf_match', self.psf_match, ['mosaic', 'group', 'none'])
        _check_option('kernel', self.kernel, ['moffat', 'gaussian'])
        _check_option('scale', self.scale, ['sources', 'calib', 'none'])
        _check_option('background', self.background, ['constant', '2d', 'none'])
        _check_option('reproject', self.reproject, ['adaptive', 'interp', 'exact'])
        for name in ['psf_margin', 'match_radius']:
            _check_positive(name, getattr(self, name))
        if self.pixel_scale is not None:
            _check_positive('pixel_scale', self.pixel_scale)
        if self.rej is not None:
            _check_positive('rej', self.rej)
