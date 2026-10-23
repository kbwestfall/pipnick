r"""
Nominal atmospheric extinction.

The extinction curve of each observatory is distributed with pipnick in
``pipnick/data/extinction``, and selected by the ``extinction`` entry of the
``[observatory]`` section of the camera configuration; see
:attr:`~pipnick.cameras.cameras.NickelCamera.extinction_file`.  The curve
gives the extinction in magnitudes per unit airmass as a function of
wavelength.

The Mt. Hamilton (Lick Observatory) curve, ``mthamextinct.dat``, is copied
unchanged from ``pypeit/data/extinction/mthamextinct.dat`` in `PypeIt
<https://github.com/pypeit/PypeIt>`__ (BSD 3-Clause; see
``licenses/README.rst``).  It covers 3200 to 10870 angstroms.  Its origin is
not documented: the PypeIt list of extinction curves notes that they were
"Generally grabbed from the Internet.  Buyer beware."

The extinction coefficient of a band (:func:`band_extinction`) is computed
with the filter response curves of the ``speclite`` package.  For now, these
are the generic Bessell UBVRI responses (``bessell-U``, etc.), not those of
the Nickel filters, so the coefficients are approximate.  They are only used
when the extinction cannot be measured, or when the user fixes it at its
nominal value; see :func:`~pipnick.pipelines.calibrate.calibrate_all`.
"""

import importlib.resources
from pathlib import Path

import numpy as np

from speclite import filters

from pipnick import log
from pipnick.pkg.exceptions import PipnickPathError


def extinction_path(name):
    """
    Return the path to an extinction curve distributed with pipnick.

    Parameters
    ----------
    name : str
        The name of the file in ``pipnick/data/extinction``.

    Returns
    -------
    :class:`pathlib.Path`
        The path to the file.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickPathError`
        Raised if the file does not exist.
    """
    path = Path(str(importlib.resources.files('pipnick') / 'data' / 'extinction' / name))
    if not path.is_file():
        raise PipnickPathError(f'No extinction curve {name} in pipnick/data/extinction.')
    return path


def load_extinction(path):
    """
    Read an extinction curve.

    Parameters
    ----------
    path : str, :class:`pathlib.Path`
        The file, with the wavelength (angstroms) and the extinction
        (magnitudes per unit airmass) in the first two columns.

    Returns
    -------
    wave, ext : :class:`numpy.ndarray`
        The wavelengths and extinction, sorted by wavelength.
    """
    wave, ext = np.loadtxt(path, usecols=(0, 1), unpack=True)
    srt = np.argsort(wave)
    return wave[srt], ext[srt]


def band_extinction(wave, ext, band):
    r"""
    Compute the extinction coefficient of a band.

    The coefficient is the extinction of a source with a flat spectrum in
    :math:`f_\lambda`, observed through the band at unit airmass:

    .. math::

        k = -2.5 \log_{10}\left[\frac{\int R(\lambda)\, \lambda\,
            10^{-0.4 k(\lambda)}\, d\lambda}{\int R(\lambda)\, \lambda\,
            d\lambda}\right],

    where :math:`R(\lambda)` is the filter response (the speclite
    ``bessell-<band>`` curve), the factor :math:`\lambda` weights by the
    number of photons, and :math:`k(\lambda)` is the extinction curve.  Where
    the band extends beyond the curve, the extinction at the nearest end of
    the curve is used, and a warning gives the fraction of the band covered.

    The effective extinction of a star depends on its spectrum, which this
    nominal value does not capture.

    Parameters
    ----------
    wave, ext : :class:`numpy.ndarray`
        The extinction curve; see :func:`load_extinction`.
    band : str
        The band: ``U``, ``B``, ``V``, ``R``, or ``I``.

    Returns
    -------
    k : float
        The extinction coefficient (magnitudes per unit airmass).
    coverage : float
        The fraction of the response-weighted band covered by the curve.
    """
    response = filters.load_filter(f'bessell-{band}')
    lam = np.asarray(response.wavelength, dtype=float)
    weight = np.asarray(response.response, dtype=float) * lam
    k_lam = np.interp(lam, wave, ext)
    k = -2.5 * np.log10(np.trapezoid(weight * 10**(-0.4 * k_lam), lam) / np.trapezoid(weight, lam))
    inside = (lam >= wave[0]) & (lam <= wave[-1])
    coverage = float(np.trapezoid(np.where(inside, weight, 0.), lam) / np.trapezoid(weight, lam))
    if coverage < 1.:
        log.warning(
            f'The extinction curve covers only {100 * coverage:.0f}% of the {band} band; the '
            'extinction at the end of the curve is used for the rest.'
        )
    return float(k), coverage


def nominal_extinction(camera, bands):
    """
    Compute the nominal extinction coefficients of the bands for a camera.

    Parameters
    ----------
    camera : :class:`~pipnick.cameras.cameras.NickelCamera`
        The camera, which sets the extinction curve of its observatory.
    bands : list
        The bands; see :func:`band_extinction`.

    Returns
    -------
    dict
        The coefficient of each band, or None if the camera has no extinction
        curve.
    """
    if camera.extinction_file is None:
        return None
    wave, ext = load_extinction(extinction_path(camera.extinction_file))
    return {band: band_extinction(wave, ext, band)[0] for band in bands}
