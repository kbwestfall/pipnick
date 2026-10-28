"""
The Moffat model of the point-spread function (PSF) of an image.

The PSF of each image is fit by the photometry (see
:func:`~pipnick.photometry.psf_photometry.build_moffat_psf_model`), and saved
in the metadata of its source table (``psf``); see
:func:`~pipnick.pipelines.photometry.photometry_all`.

The model is either circular or elliptical.  Its FWHM (:attr:`PSFModel.fwhm`)
is the **mean** of the FWHM along its major and minor axes, which are also
provided (:attr:`PSFModel.fwhm_major`, :attr:`PSFModel.fwhm_minor`); the two
are equal for a circular model.  Steps that need the extent of the PSF in
every direction, e.g., the PSF matching of the mosaics, use the major axis.
"""

import numpy as np

from photutils.psf import MoffatPSF

from pipnick.photometry.moffat import (
    EllipticalMoffatPSF, FitEllipticalMoffat2D, FitMoffat2D, fwhm_to_gamma, gamma_to_fwhm
)
from pipnick.pkg.exceptions import PipnickParameterError


class PSFModel:
    """
    A circular or elliptical Moffat model of a PSF.

    The parameters are those of the fit to the PSF; see
    :class:`~pipnick.photometry.moffat.FitMoffat2D` and
    :class:`~pipnick.photometry.moffat.FitEllipticalMoffat2D`.  The model is
    defined by its shape; its position, amplitude, and background are not
    used.

    Parameters
    ----------
    fittype : str
        ``'circ'`` or ``'ellip'``.
    par : array-like
        The parameters: ``x0``, ``y0``, amplitude, ``gamma``, ``alpha``, and
        background for a circular model, or ``x0``, ``y0``, amplitude,
        ``gamma1``, ``gamma2``, ``phi``, ``alpha``, and background for an
        elliptical model.  The core widths (``gamma``) are in pixels, ``phi``
        is in radians, and ``alpha`` is the power index.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if the fit type is unknown, or the number of parameters does
        not match it.
    """
    def __init__(self, fittype, par):
        npar = {'circ': 6, 'ellip': 8}
        if fittype not in npar:
            raise PipnickParameterError(
                f'Unknown PSF fit type {fittype}; must be circ or ellip.'
            )
        if len(par) != npar[fittype]:
            raise PipnickParameterError(
                f'A {fittype} PSF model has {npar[fittype]} parameters, not {len(par)}.'
            )
        self.fittype = fittype
        self.par = tuple(float(p) for p in par)

    @classmethod
    def from_fit(cls, fit):
        """
        Construct the model from a fit to a PSF.

        Parameters
        ----------
        fit : object
            The fit, a :class:`~pipnick.photometry.moffat.FitMoffat2D` or
            :class:`~pipnick.photometry.moffat.FitEllipticalMoffat2D` object.

        Returns
        -------
        :class:`PSFModel`
            The model.
        """
        if isinstance(fit, FitEllipticalMoffat2D):
            return cls('ellip', fit.par)
        if isinstance(fit, FitMoffat2D):
            return cls('circ', fit.par)
        raise PipnickParameterError(f'Cannot construct a PSF model from {type(fit).__name__}.')

    @classmethod
    def from_meta(cls, meta):
        """
        Construct the model from the metadata of a source table.

        Parameters
        ----------
        meta : dict
            The ``psf`` entry of the metadata of a source table, with the
            ``fittype`` and ``par``; see :meth:`to_meta`.

        Returns
        -------
        :class:`PSFModel`
            The model.
        """
        return cls(meta['fittype'], meta['par'])

    @classmethod
    def circular(cls, fwhm, beta):
        """
        Construct a circular model with a given FWHM and power index.

        Parameters
        ----------
        fwhm : float
            The FWHM in pixels.
        beta : float
            The power index.

        Returns
        -------
        :class:`PSFModel`
            The model, centered at (0, 0) with unit amplitude.
        """
        return cls('circ', [0., 0., 1., fwhm_to_gamma(fwhm, beta), beta, 0.])

    @property
    def beta(self):
        """The power index."""
        return self.par[4] if self.fittype == 'circ' else self.par[6]

    @property
    def gamma(self):
        """The core widths (pixels) along the two axes of the model."""
        return (self.par[3], self.par[3]) if self.fittype == 'circ' else self.par[3:5]

    @property
    def phi(self):
        """The rotation of the model (radians; 0 for a circular model)."""
        return 0. if self.fittype == 'circ' else self.par[5]

    @property
    def fwhm_major(self):
        """The FWHM (pixels) along the major axis."""
        return float(gamma_to_fwhm(max(np.abs(self.gamma)), self.beta))

    @property
    def fwhm_minor(self):
        """The FWHM (pixels) along the minor axis."""
        return float(gamma_to_fwhm(min(np.abs(self.gamma)), self.beta))

    @property
    def fwhm(self):
        """The FWHM (pixels): the mean of the FWHM along the two axes."""
        return float(gamma_to_fwhm(np.mean(np.abs(self.gamma)), self.beta))

    def to_photutils(self, **kwargs):
        """
        Return the ``photutils`` PSF model with this shape.

        Parameters
        ----------
        **kwargs
            Other parameters of the model, e.g., ``flux``, ``x_0``, and
            ``y_0``.

        Returns
        -------
        :class:`astropy.modeling.Fittable2DModel`
            The model: a :class:`photutils.psf.MoffatPSF` (which calls the core
            width ``alpha`` and the power index ``beta``) or an
            :class:`~pipnick.photometry.moffat.EllipticalMoffatPSF`.
        """
        if self.fittype == 'circ':
            return MoffatPSF(alpha=self.par[3], beta=self.par[4], **kwargs)
        return EllipticalMoffatPSF(
            gamma1=self.par[3], gamma2=self.par[4], phi=self.par[5], alpha=self.par[6], **kwargs
        )

    def render(self, size, scale=1.):
        """
        Render the model, normalized to unit sum, at the center of an image.

        Parameters
        ----------
        size : int
            The size of the (square) image; must be odd.
        scale : float, optional
            The ratio of the pixel size of the model to that of the image,
            used to render a model fit to an image with a different pixel
            scale.

        Returns
        -------
        :class:`numpy.ndarray`
            The PSF image.
        """
        yy, xx = np.mgrid[:size, :size]
        c = (size - 1) / 2
        if self.fittype == 'circ':
            model = MoffatPSF(flux=1., x_0=c, y_0=c, alpha=self.par[3] * scale, beta=self.par[4])
        else:
            model = EllipticalMoffatPSF(
                flux=1., x_0=c, y_0=c, gamma1=self.par[3] * scale, gamma2=self.par[4] * scale,
                phi=self.par[5], alpha=self.par[6]
            )
        psf = model(xx, yy)
        return psf / np.sum(psf)

    def to_meta(self, pixel_scale, nstack):
        """
        Return the metadata saved with a source table.

        Parameters
        ----------
        pixel_scale : float
            The pixel scale of the image (arcsec per pixel).
        nstack : int
            The number of sources used to fit the model.

        Returns
        -------
        dict
            The fit type (``fittype``), the parameters (``par``), the FWHM in
            pixels and arcsec (``fwhm_pix``, ``fwhm_arcsec``; see
            :attr:`fwhm`), and ``nstack``.
        """
        return {
            'fittype': self.fittype, 'par': list(self.par), 'fwhm_pix': self.fwhm,
            'fwhm_arcsec': float(self.fwhm * pixel_scale), 'nstack': int(nstack)
        }
