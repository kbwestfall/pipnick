r"""
Matching the point-spread functions (PSFs) of images.

Images to be combined are matched to a common PSF by convolving each with a
kernel :math:`K` such that :math:`P \ast K = T`, where :math:`P` is the PSF of
the image and :math:`T` is the target PSF.  Convolution can only broaden a
PSF, so the target must be broader than every PSF, and its wings must be at
least as extended.  The target is a circular Moffat function with a FWHM
slightly larger than the largest FWHM of the images.  For an elliptical PSF,
the FWHM along its major axis is used, so that the kernel broadens the PSF in
every direction (see :attr:`~pipnick.datamodel.psf.PSFModel.fwhm_major`).
Its power index is the median of those of the images: in principle, the
smallest power index (the most extended wings) should be used, but the power
indices of the PSF models are poorly constrained (they are measured from
stacks of a few stars, and are degenerate with the FWHM), and the smallest
selects the noisiest measurement.

The PSF of each image is the Moffat model fit by the photometry (see
:func:`~pipnick.photometry.psf_photometry.build_moffat_psf_model` and
:class:`~pipnick.datamodel.psf.PSFModel`).  The
kernel is computed with :func:`photutils.psf.matching.make_kernel` (see
:func:`matching_kernel`), or, if that fails, approximated by a Gaussian (see
:func:`gaussian_kernel`).

The images are convolved with *normalized convolution* (see
:func:`normalized_convolve`), so that masked pixels do not contribute.  The
variance is convolved with the same kernel, which gives the variance of sums
over apertures much larger than the kernel, but overestimates the variance of
single pixels, whose noise becomes correlated; see
:func:`normalized_convolve`.
"""

import numpy as np

from astropy.convolution import Gaussian2DKernel
from photutils.psf_matching import SplitCosineBellWindow, make_kernel
from scipy.signal import fftconvolve

from pipnick import log
from pipnick.datamodel.psf import PSFModel


def kernel_size(fwhm, factor=6.):
    """
    Return the (odd) size of the images used for the PSFs and kernels:
    ``factor`` times the FWHM (pixels), and at least 11.
    """
    size = max(11, int(np.ceil(factor * fwhm)))
    return size + 1 - size % 2


def gaussian_kernel(fwhm, target_fwhm, size):
    r"""
    Return the Gaussian kernel that broadens a Gaussian PSF with FWHM
    ``fwhm`` to one with FWHM ``target_fwhm``.

    The kernel has :math:`\sigma_k^2 = \sigma_t^2 - \sigma^2`.  If the target
    is not broader, the kernel is a delta function.

    Parameters
    ----------
    fwhm, target_fwhm : float
        The FWHM of the PSF and of the target, in pixels.
    size : int
        The size of the kernel; must be odd.

    Returns
    -------
    :class:`numpy.ndarray`
        The kernel, normalized to unit sum.
    """
    var = (target_fwhm**2 - fwhm**2) / (8 * np.log(2))
    if var <= 0:
        kernel = np.zeros((size, size), dtype=float)
        kernel[size // 2, size // 2] = 1.
        return kernel
    kernel = Gaussian2DKernel(np.sqrt(var), x_size=size, y_size=size).array
    return kernel / np.sum(kernel)


def matching_kernel(psf, target, max_residual=0.05):
    r"""
    Compute the kernel that matches a PSF to a target PSF.

    The kernel is computed by :func:`photutils.psf.matching.make_kernel` with a
    :class:`photutils.psf.matching.SplitCosineBellWindow`, and normalized to
    unit sum.  It is accepted only if the convolved PSF matches the target to
    within ``max_residual`` (the root-sum-square of the difference, relative
    to that of the target), and it does not amplify the noise, i.e.,
    :math:`\sum K^2 < 1`.  The kernel of a PSF that is almost as broad as the
    target is close to a delta function, with negative ringing; this is
    accurate, and reduces the noise, so the negative values are not a reason
    to reject it.

    Parameters
    ----------
    psf, target : :class:`numpy.ndarray`
        The PSF and target images, with the same odd size, normalized to unit
        sum.
    max_residual : float, optional
        The largest acceptable relative residual.

    Returns
    -------
    kernel : :class:`numpy.ndarray`
        The kernel, or None if it is not acceptable.
    residual : float
        The relative residual of the convolved PSF.
    """
    kernel = make_kernel(psf, target, window=SplitCosineBellWindow(alpha=0.4, beta=0.3))
    total = np.sum(kernel)
    if total <= 0:
        return None, 1.
    kernel = kernel / total
    residual = float(
        np.sqrt(np.sum((fftconvolve(psf, kernel, mode='same') - target)**2))
        / np.sqrt(np.sum(target**2))
    )
    if residual > max_residual or np.sum(kernel**2) >= 1.:
        return None, residual
    return kernel, residual


def circular_taper(kernel, start=0.85):
    """
    Taper a kernel to zero outside a circle, and renormalize it.

    A matching kernel computed on a square image extends to its edges, which
    gives a square halo around bright stars.  The kernel is multiplied by a
    window that is 1 within ``start`` times the radius of the inscribed
    circle, falls to 0 at the radius as a cosine, and is 0 outside it.

    Parameters
    ----------
    kernel : :class:`numpy.ndarray`
        The (square, odd-sized) kernel.
    start : float, optional
        The fraction of the radius where the taper starts.

    Returns
    -------
    :class:`numpy.ndarray`
        The tapered kernel, normalized to unit sum.
    """
    size = kernel.shape[0]
    c = (size - 1) / 2
    yy, xx = np.mgrid[:size, :size]
    r = np.hypot(xx - c, yy - c) / c
    x = np.clip((r - start) / (1 - start), 0., 1.)
    window = 0.5 * (1 + np.cos(np.pi * x))
    tapered = kernel * window
    return tapered / np.sum(tapered)


def image_kernel(psf, target_fwhm, target_beta, method='moffat'):
    """
    Compute the kernel that matches the PSF of an image to the target.

    Parameters
    ----------
    psf : :class:`~pipnick.datamodel.psf.PSFModel`
        The PSF model of the image.
    target_fwhm : float
        The FWHM of the target, in pixels of the image.
    target_beta : float
        The power index of the target.
    method : str, optional
        ``'moffat'``: the kernel matching the Moffat models (see
        :func:`matching_kernel`), with a Gaussian kernel as the fallback if it
        fails; or ``'gaussian'``: the Gaussian kernel (see
        :func:`gaussian_kernel`).

    Returns
    -------
    kernel : :class:`numpy.ndarray`
        The kernel.  A Moffat kernel is tapered to a circle; see
        :func:`circular_taper`.
    used : str
        The kernel used: ``'moffat'``, ``'gaussian'``, or ``'gaussian
        (fallback)'``.
    ee_ratio : float
        The accuracy of the matching: the ratio of the encircled energy of the
        matched PSF to that of the target within 1.5 times the target FWHM;
        see :func:`encircled_energy_ratio`.
    """
    # The kernel is tapered to a circle (see circular_taper), so its size is
    # increased to keep the extent of its wings
    size = kernel_size(target_fwhm, factor=10.)
    # The extent of the PSF in every direction; see PSFModel
    fwhm = psf.fwhm_major
    image = psf.render(size)
    target = PSFModel.circular(target_fwhm, target_beta).render(size)
    kernel = None
    if method == 'moffat':
        kernel, residual = matching_kernel(image, target)
        if kernel is None:
            log.warning(
                f'The Moffat matching kernel failed (residual {residual:.3f}); using a Gaussian '
                'kernel instead.'
            )
            used = 'gaussian (fallback)'
        else:
            kernel = circular_taper(kernel)
            used = 'moffat'
    else:
        used = 'gaussian'
    if kernel is None:
        kernel = gaussian_kernel(fwhm, target_fwhm, size)
    return kernel, used, encircled_energy_ratio(image, kernel, target, 1.5 * target_fwhm)


def encircled_energy_ratio(psf, kernel, target, radius):
    """
    Return the ratio of the encircled energy of a PSF convolved with a kernel
    to that of the target, within ``radius`` (pixels) of the center: a
    measure of the accuracy of the matching, for aperture photometry.
    """
    size = psf.shape[0]
    c = (size - 1) / 2
    yy, xx = np.mgrid[:size, :size]
    inside = np.hypot(xx - c, yy - c) < radius
    return float(np.sum(fftconvolve(psf, kernel, mode='same')[inside]) / np.sum(target[inside]))


def normalized_convolve(data, var, good, kernel, min_weight=0.9):
    r"""
    Convolve an image, its variance, and its good-pixel map with a kernel,
    ignoring the masked pixels.

    The convolved image is

    .. math::

        d' = \frac{(d\, g) \ast K}{g \ast K},

    where :math:`g` is 1 for good pixels and 0 for masked pixels, so that the
    masked pixels do not contribute.  Pixels with :math:`g \ast K` less than
    ``min_weight`` (i.e., where too much of the kernel falls on masked pixels)
    are masked.

    The variance is convolved in the same way, with :math:`K` rather than
    :math:`K^2`.  The result is a *variance density*: summing it over an
    aperture much larger than the kernel gives the variance of the sum of the
    convolved image over the aperture, but it overestimates the variance of a
    single pixel, because the noise of neighboring pixels becomes correlated.

    Parameters
    ----------
    data, var : :class:`numpy.ndarray`
        The image and its variance.
    good : :class:`numpy.ndarray`
        Boolean array selecting the good pixels.
    kernel : :class:`numpy.ndarray`
        The kernel, normalized to unit sum.
    min_weight : float, optional
        The smallest fraction of the kernel that must fall on good pixels.

    Returns
    -------
    data, var : :class:`numpy.ndarray`
        The convolved image and variance, 0 at the masked pixels.
    good : :class:`numpy.ndarray`
        Boolean array selecting the good pixels of the convolved image.
    """
    g = np.asarray(good, dtype=float)
    weight = fftconvolve(g, kernel, mode='same')
    _good = np.asarray(good, dtype=bool) & (weight > min_weight)
    _weight = np.where(_good, weight, 1.)
    _data = fftconvolve(np.where(good, data, 0.), kernel, mode='same') / _weight
    _var = fftconvolve(np.where(good, var, 0.), kernel, mode='same') / _weight
    return np.where(_good, _data, 0.), np.where(_good, np.maximum(_var, 0.), 0.), _good
