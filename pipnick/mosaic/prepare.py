"""
Prepare the images of a mosaic, and measure their flux scales.

Each image is converted to electrons per second, with its background
subtracted (see :func:`prepare_image`), and the magnitude offsets of the
images are measured from the sources common to overlapping images, or from
their photometric calibrations (see :func:`measure_offsets` and
:mod:`pipnick.mosaic.scaling`).
"""

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from photutils.aperture import ApertureStats, CircularAnnulus
from photutils.segmentation import detect_sources, detect_threshold

from pipnick import log
from pipnick.image.io import invalid_error, read_reduced
from pipnick.mosaic.scaling import solve_offsets
from pipnick.photometry.aperture_photometry import aperture_analysis
from pipnick.photometry.psf_photometry import subtract_background
from pipnick.pkg.exceptions import PipnickDataModelError
from pipnick.utils.crossmatch import match_sky


def _background(data, good, method):
    """
    Return the background of an image: ``'constant'`` (the sigma-clipped
    median of the pixels without sources, which are found by image
    segmentation and dilated), ``'2d'`` (see
    :func:`~pipnick.photometry.psf_photometry.subtract_background`, with
    256-pixel boxes), or ``'none'`` (0).
    """
    if method == 'none':
        return 0.
    mask = np.logical_not(good)
    if method == '2d':
        return data - subtract_background(data, mask=mask, box_size=256)
    threshold = detect_threshold(data, 2., mask=mask)
    segm = detect_sources(data, threshold, 10, mask=mask)
    sources = np.zeros(data.shape, dtype=bool) if segm is None else segm.make_source_mask(size=11)
    _, median, _ = sigma_clipped_stats(data, mask=mask | sources, sigma=3.)
    return median



def prepare_image(image, background='constant'):
    """
    Read an image, convert it to electrons per second, and subtract its
    background.

    Returns
    -------
    data, var : :class:`numpy.ndarray`
        The image and its variance, 0 at the masked pixels.
    good : :class:`numpy.ndarray`
        Boolean array selecting the good pixels.
    """
    data, mask, err, _ = read_reduced(image.path)
    if err is None:
        log.warning(f'{image.name} has no uncertainty; using unit errors.')
        err = np.ones(data.shape, dtype=float)
    good = np.logical_not(mask | invalid_error(err)) & np.isfinite(data)
    t = image.exptime
    data = np.where(good, data, 0.) / t
    var = np.where(good, err, 0.)**2 / t**2
    data = np.where(good, data - _background(data, good, background), 0.)
    return data, var, good



def _source_mags(image, data, var, good, target_pix, settings):
    """
    Return the positions and instrumental magnitudes of the sources of an
    image used to measure the flux scale.

    The sources are those of the photometry table that are not flagged, have
    an unmasked aperture flux (so no saturated or masked pixels), and have
    celestial coordinates.  With PSF matching, their fluxes are measured on
    the matched image in apertures of radius 1.5 times the target FWHM, with
    the local background (the sigma-clipped median in an annulus from 3 to 5
    times the target FWHM) subtracted, so that errors in the background of the
    image do not bias the scale; otherwise, their PSF fluxes are used.
    Sources with a signal-to-noise ratio below ``settings['min_snr']``, or
    without a local background, are excluded.
    """
    phot = image.photometry(require=False)
    ok = phot.usable()
    x = np.asarray(phot['x'], dtype=float)[ok]
    y = np.asarray(phot['y'], dtype=float)[ok]
    ra = np.ma.filled(np.ma.asarray(phot['ra'], dtype=float), 0.)[ok]
    dec = np.ma.filled(np.ma.asarray(phot['dec'], dtype=float), 0.)[ok]
    if settings['psf_match'] != 'none':
        mask = np.logical_not(good)
        if len(x) > 0:
            annulus = CircularAnnulus(np.column_stack([x, y]), 3 * target_pix, 5 * target_pix)
            # NOTE: The median is NaN for an annulus without good pixels; these
            # sources are flagged and excluded immediately
            bkg = np.atleast_1d(np.asarray(ApertureStats(data, annulus, mask=mask).median))
        else:
            bkg = np.zeros(0)
        has_bkg = np.isfinite(bkg)
        flux, err = aperture_analysis(
            data, x, y, 1.5 * target_pix, local_bkg=np.where(has_bkg, bkg, 0.),
            error=np.sqrt(var), mask=mask
        )
        flux = np.where(has_bkg, np.ma.filled(flux, 0.), 0.)
        err = np.where(has_bkg, np.ma.filled(err, 0.), 0.)
    else:
        t = image.exptime
        flux = np.ma.filled(np.ma.asarray(phot['flux_psf'], dtype=float), 0.)[ok] / t
        err = np.ma.filled(np.ma.asarray(phot['flux_psf_err'], dtype=float), 0.)[ok] / t
    keep = (flux > 0) & (err > 0) & (flux > settings['min_snr'] * err)
    _flux = np.where(keep, flux, 1.)
    mag = -2.5 * np.log10(_flux)
    mag_err = 2.5 / np.log(10.) * np.where(keep, err, 0.) / _flux
    return ra[keep], dec[keep], mag[keep], mag_err[keep]



def _reference(images, reference):
    """
    Return the index of the reference image: the image named ``reference``
    (file name or stem), if it is in the mosaic, or the image with the lowest
    airmass.
    """
    if reference is not None:
        for i, image in enumerate(images):
            if reference in [image.name, image.path.stem]:
                return i
    airmass = [np.inf if image.airmass is None else image.airmass for image in images]
    return int(np.argmin(airmass))



def measure_offsets(images, prepared, targets, settings):
    """
    Measure the magnitude offsets of the images of a mosaic with respect to
    the reference image.

    Parameters
    ----------
    images : list
        The images of the mosaic; see :func:`~pipnick.pipelines.mosaic.collect_images`.
    prepared : list
        The matched images: tuples of the image, variance, and good-pixel map.
    targets : list
        The target FWHM of each image, in its pixels.
    settings : dict
        The mosaic settings; see :func:`~pipnick.pipelines.mosaic.mosaic_all`.

    Returns
    -------
    offset, error : :class:`numpy.ndarray`
        The magnitude offset of each image (it is multiplied by
        :math:`10^{0.4 z}`) and its error.
    scaled : :class:`numpy.ndarray`
        Flags the images whose offset was measured.
    nused : :class:`numpy.ndarray`
        The number of magnitude differences used for each image.
    reference : int
        The index of the reference image.
    """
    n = len(images)
    reference = _reference(images, settings['reference'])
    if settings['scale'] == 'none' or n == 1:
        return np.zeros(n), np.zeros(n), np.ones(n, dtype=bool), np.zeros(n, dtype=int), reference
    if settings['scale'] == 'calib':
        zp = [image.header.get('PNCALZP') for image in images]
        if any(z is None for z in zp):
            raise PipnickDataModelError(
                'Not all images have a photometric calibration (PNCALZP); run pipnick_calibrate '
                '--apply image, or use --scale sources.'
            )
        offset = np.array([zp[reference] - z for z in zp], dtype=float)
        return offset, np.zeros(n), np.ones(n, dtype=bool), np.zeros(n, dtype=int), reference
    sources = [
        _source_mags(image, d, v, g, t, settings)
        for image, (d, v, g), t in zip(images, prepared, targets)
    ]
    measurements = []
    for i in range(n):
        for j in range(i + 1, n):
            i1, i2, _ = match_sky(
                SkyCoord(sources[i][0], sources[i][1], unit='deg'),
                SkyCoord(sources[j][0], sources[j][1], unit='deg'), settings['match_radius']
            )
            if len(i1) == 0:
                continue
            dm = sources[i][2][i1] - sources[j][2][i2]
            err = np.sqrt(sources[i][3][i1]**2 + sources[j][3][i2]**2)
            measurements.append((i, j, dm, err))
    offset, error, scaled, nused = solve_offsets(
        n, measurements, reference=reference, min_common=settings['min_common']
    )
    for image, s in zip(images, scaled):
        if not s:
            log.warning(
                f'{image.name} shares too few sources with the other images (fewer than '
                f'{settings["min_common"]} per pair) to be scaled; it is not scaled.'
            )
    return offset, error, scaled, nused, reference

