r"""
Mosaics of the reduced science images.

The images to combine are selected by the ``mosaic`` and ``wcs_group``
columns of the per-setup metadata tables (see
:data:`~pipnick.metadata.mosaic_columns`): images with the same ``mosaic``
label are combined, and all of the mosaics with the same ``wcs_group`` label
share an identical output WCS and shape (e.g., for colour images).  The
labels are global, so a mosaic can combine the images of several tables,
e.g., of different setups.

A mosaic can have a single image, which is then only reprojected onto the WCS
of its group, e.g., so that the images of an object in several filters are
sampled on the same pixel grid for a colour image (see
:mod:`pipnick.utils.rgb`).

The images must have astrometric solutions (``pipnick_astrometry``) and,
unless both the PSF matching and the flux scaling are off, current
photometry source tables (``pipnick_photometry``).  Each mosaic is
constructed by:

    1. converting each image to electrons per second, and subtracting its
       background;
    2. convolving each image to a common PSF (see :mod:`pipnick.mosaic.psfmatch`);
    3. scaling each image to the flux scale of a reference image, using the
       sources common to overlapping images (see :mod:`pipnick.mosaic.scaling`);
    4. reprojecting the images onto the output WCS of the WCS group and
       combining them, rejecting outliers (see :mod:`pipnick.mosaic.coadd`).

For now, a mosaic must combine images from a single night; see
:func:`mosaic_all`.
"""

from pathlib import Path

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.visualization import ZScaleInterval
from matplotlib.figure import Figure
from photutils.aperture import ApertureStats, CircularAnnulus
from photutils.segmentation import detect_sources, detect_threshold

from pipnick import log
from pipnick.astrometry.fit import wcs_from_header
from pipnick.datamodel.dataset import DataSet
from pipnick.datamodel.psf import PSFModel
from pipnick.image.io import invalid_error, read_reduced
from pipnick.mosaic.coadd import (
    combine, footprint_bbox, output_wcs, overlap_differences, reproject_image
)
from pipnick.mosaic.psfmatch import image_kernel, normalized_convolve
from pipnick.mosaic.scaling import solve_offsets
from pipnick.photometry.aperture_photometry import aperture_analysis
from pipnick.photometry.psf_photometry import build_moffat_psf_model, subtract_background
from pipnick.pkg.exceptions import (
    PipnickDataModelError, PipnickParameterError, PipnickPathError
)
from pipnick.pkg.paths import step_dir
from pipnick.pkg.provenance import image_hash, json_hash, solution_hash
from pipnick.pkg.report import StepReport
from pipnick.pkg.version import version
from pipnick.utils.crossmatch import match_sky


def collect_images(tables, require_photometry=True, report=None):
    """
    Collect the reduced science images assigned to mosaics.

    Parameters
    ----------
    tables : list
        The per-setup metadata files.
    require_photometry : bool, optional
        Require a current photometry source table for each image.
    report : :class:`~pipnick.pkg.report.StepReport`, optional
        Records the pipnick version of each table.

    Returns
    -------
    list
        The :class:`~pipnick.datamodel.image.ReducedImage` objects with a
        mosaic label; see :func:`mosaic_label` and :func:`wcs_group`.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the labels are invalid (see
        :func:`~pipnick.datamodel.dataset.check_mosaic_groups`), an image is
        not from the setup of its table or has no astrometric solution,
        exposure time, or MJD, or its photometry table is missing (when
        required) or out of date.
    """
    images = []
    for table in tables:
        dataset = DataSet.from_file(table, check_mosaics=True)
        if report is not None:
            report.check_version(dataset.metafile.name, dataset.table.meta.get('pipnick_version'))
        if 'mosaic' not in dataset.table.colnames or np.all(dataset.table['mosaic'] == ''):
            log.warning(f'{dataset.metafile.name} defines no mosaics.')
            continue
        for image in dataset.science_images():
            if len(mosaic_label(image)) == 0:
                continue
            if not image.has_solution:
                raise PipnickDataModelError(
                    f'{image.name} has no astrometric solution; run pipnick_astrometry.'
                )
            # Check that the photometry table exists (if required) and is up to
            # date
            image.photometry(require=require_photometry)
            if image.exptime is None or image.mjd is None:
                raise PipnickDataModelError(f'{image.name} has no exposure time or MJD.')
            images.append(image)
    return images


def mosaic_label(image):
    """
    Return the mosaic label of an image; see
    :meth:`~pipnick.datamodel.dataset.DataSet.mosaic_labels`.
    """
    return image.dataset.mosaic_labels(image.product)[0]


def wcs_group(image):
    """
    Return the WCS group of an image, which defaults to its mosaic label; see
    :meth:`~pipnick.datamodel.dataset.DataSet.mosaic_labels`.
    """
    return image.dataset.mosaic_labels(image.product)[1]


def _validate(images):
    """
    Check the mosaics across all of the tables: each mosaic has one WCS group
    and one night.  Mixed filters are reported with a warning.

    Returns
    -------
    dict
        The images of each mosaic.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if a mosaic has several WCS groups or nights.
    """
    mosaics = {}
    for image in images:
        mosaics.setdefault(mosaic_label(image), []).append(image)
    for label, members in mosaics.items():
        groups = sorted({wcs_group(m) for m in members})
        if len(groups) > 1:
            raise PipnickDataModelError(
                f'Mosaic {label} has several WCS groups in different tables: '
                f'{", ".join(groups)}.'
            )
        nights = sorted({m.night for m in members})
        if len(nights) > 1:
            raise PipnickDataModelError(
                f'Mosaic {label} combines images from several nights ({", ".join(nights)}); '
                'mosaics of several nights are not yet supported.'
            )
        filters = sorted({str(m.filter) for m in members})
        if len(filters) > 1:
            log.warning(
                f'Mosaic {label} combines images with different filters: {", ".join(filters)}.'
            )
        setups = sorted({m.dataset.tag for m in members})
        log.info(
            f'Mosaic {label}: {len(members)} images; WCS group {groups[0]}; night {nights[0]}; '
            f'setups {", ".join(setups)}.'
        )
    return mosaics


def _psf(image):
    """
    Return the PSF model of an image, from its photometry table, with its FWHM
    along the major axis (arcsec; see
    :attr:`~pipnick.datamodel.psf.PSFModel.fwhm_major`) and its power index.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the image has no photometry table or PSF model.
    """
    phot = image.photometry(require=False)
    if phot is None or phot.meta.get('psf', {}).get('fwhm_pix', -1) <= 0:
        raise PipnickDataModelError(
            f'{image.name} has no PSF model; run pipnick_photometry, or turn off the '
            'PSF matching.'
        )
    psf = PSFModel.from_meta(phot.meta['psf'])
    return psf, psf.fwhm_major * image.pixel_scale, psf.beta


def psf_target(images, margin=1.05, outlier=1.5):
    """
    Return the target PSF of a set of images: the FWHM (arcsec), ``margin``
    times the largest FWHM, and the median power index; see
    :mod:`pipnick.mosaic.psfmatch`.

    A warning is issued for images whose FWHM is more than ``outlier`` times
    the median, because they degrade the resolution of all of the others.
    """
    psfs = [_psf(image) for image in images]
    fwhm = np.array([p[1] for p in psfs])
    median = np.median(fwhm)
    for image, f in zip(images, fwhm):
        if f > outlier * median:
            log.warning(
                f'The PSF FWHM of {image.name} ({f:.2f} arcsec) is much larger than the '
                f'median ({median:.2f} arcsec), and sets the target of the PSF matching.  '
                'Consider excluding it from the mosaic.'
            )
    return margin * float(np.max(fwhm)), float(np.median([p[2] for p in psfs]))


def group_wcs(images, wcs_file, pixel_scale=None, overwrite=False):
    """
    Return the output WCS and shape of a WCS group.

    The WCS covers all of the images in the group (see
    :func:`~pipnick.mosaic.coadd.output_wcs`).  It is saved to ``wcs_file``
    (a header), with the hash of the group membership (``PNGRPHSH``), and
    reused while the membership is unchanged, so that every mosaic of the
    group has the identical WCS and shape.

    Returns
    -------
    wcs : :class:`astropy.wcs.WCS`
        The WCS.
    shape : tuple
        The shape.
    ghash : str
        The hash of the group membership and pixel scale.
    """
    ghash = json_hash(
        {
            'images': sorted(
                [str(image.path), solution_hash(image.header)] for image in images
            ),
            'pixel_scale': pixel_scale
        }
    )
    if wcs_file.is_file() and not overwrite:
        hdr = fits.getheader(wcs_file)
        if hdr.get('PNGRPHSH') == ghash:
            return wcs_from_header(hdr), (int(hdr['PNNAXIS2']), int(hdr['PNNAXIS1'])), ghash
        log.info(f'The membership of WCS group {wcs_file.stem} changed; recomputing its WCS.')
    wcs, shape = output_wcs(
        [image.wcs for image in images], [image.shape for image in images],
        pixel_scale=pixel_scale
    )
    hdr = wcs.to_header()
    hdr['PNNAXIS1'] = (shape[1], 'Number of columns of the mosaics')
    hdr['PNNAXIS2'] = (shape[0], 'Number of rows of the mosaics')
    hdr['PNGRPHSH'] = (ghash, 'Hash of the WCS group membership')
    hdr['PNVERS'] = (version, 'pipnick version (mosaic)')
    wcs_file.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(header=hdr).writeto(wcs_file, overwrite=True)
    log.info(f'Wrote the WCS of group {wcs_file.stem} ({shape[1]} x {shape[0]} pixels).')
    return wcs, shape, ghash


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
    ok = (
        (np.asarray(phot['flags']) == 0)
        & np.logical_not(np.ma.getmaskarray(np.ma.asarray(phot['flux_aper'])))
        & np.logical_not(np.ma.getmaskarray(np.ma.asarray(phot['ra'])))
    )
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


def measure_offsets(images, prepared, targets, settings):
    """
    Measure the magnitude offsets of the images of a mosaic with respect to
    the reference image.

    Parameters
    ----------
    images : list
        The images of the mosaic; see :func:`collect_images`.
    prepared : list
        The matched images: tuples of the image, variance, and good-pixel map.
    targets : list
        The target FWHM of each image, in its pixels.
    settings : dict
        The mosaic settings; see :func:`mosaic_all`.

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


def build_mosaic(images, wcs_out, shape_out, target, settings):
    """
    Construct a mosaic.

    Parameters
    ----------
    images : list
        The images of the mosaic; see :func:`collect_images`.
    wcs_out : :class:`astropy.wcs.WCS`
        The output WCS.
    shape_out : tuple
        The output shape.
    target : tuple
        The target FWHM (arcsec) and power index of the PSF matching, or None
        if the PSFs are not matched.
    settings : dict
        The mosaic settings; see :func:`mosaic_all`.

    Returns
    -------
    hdus : :class:`astropy.io.fits.HDUList`
        The mosaic; see :func:`mosaic_all`.
    inputs : :class:`astropy.table.Table`
        The table of the input images.
    """
    prepared = []
    targets = []
    kernels = []
    ee_ratios = []
    for image in images:
        data, var, good = prepare_image(image, background=settings['background'])
        if target is not None:
            psf, _, _ = _psf(image)
            target_pix = target[0] / image.pixel_scale
            kernel, used, ee_ratio = image_kernel(
                psf, target_pix, target[1], method=settings['kernel'], name=image.name
            )
            data, var, good = normalized_convolve(data, var, good, kernel)
            kernels.append(used)
            ee_ratios.append(ee_ratio)
        else:
            # The FWHM (pixels) used for the scaling apertures, which are not
            # used without PSF matching (the PSF fluxes are used instead)
            target_pix = 0.
            kernels.append('none')
            ee_ratios.append(-1.)
        prepared.append((data, var, good))
        targets.append(target_pix)

    offset, offset_err, scaled, nused, ref = measure_offsets(images, prepared, targets, settings)
    weights = []
    reprojected = []
    achieved = []
    for image, (data, var, good), z, t in zip(images, prepared, offset, targets):
        s = 10**(0.4 * z)
        data = data * s
        var = var * s**2
        achieved.append(
            -1. if target is None else _achieved_fwhm(image, data, good, t, target[1])
        )
        weights.append(1. / float(np.median(var[good])))
        bbox = footprint_bbox(image.wcs, image.shape, wcs_out, shape_out)
        if bbox is None:
            raise PipnickDataModelError(f'{image.name} does not overlap the mosaic.')
        sub = wcs_out[bbox]
        shape = (bbox[0].stop - bbox[0].start, bbox[1].stop - bbox[1].start)
        d, v, g = reproject_image(
            data, var, good, image.wcs, sub, shape, method=settings['reproject']
        )
        reprojected.append((bbox, d, v, g))
    _check_overlap(images, reprojected)
    levels = np.zeros(len(images), dtype=float)
    if settings['match_background'] and len(images) > 1:
        levels, _, matched, _ = solve_offsets(
            len(images), overlap_differences(reprojected), reference=ref, floor=0.,
            min_common=1
        )
        reprojected = [
            (bbox, np.where(g, d - b, 0.), v, g)
            for (bbox, d, v, g), b in zip(reprojected, levels)
        ]
        log.info(
            'Matched the background levels (e-/s): '
            + ', '.join(f'{b:+.4f}' for b in levels)
        )

    mean, var, coverage, used, nrej = combine(
        reprojected, shape_out, weights, rej=settings['rej'], rej_frac=settings['rej_frac']
    )
    expmap = np.zeros(shape_out, dtype=float)
    nreject = np.zeros(shape_out, dtype=int)
    for (bbox, _, _, g), u, image in zip(reprojected, used, images):
        expmap[bbox] += u * image.exptime
        nreject[bbox] += g & np.logical_not(u)
    if nrej > 0:
        log.info(f'Rejected {nrej} outlying pixels.')

    inputs = Table()
    inputs['image'] = [image.name for image in images]
    inputs['table'] = [image.dataset.metafile.name for image in images]
    inputs['setup'] = [image.dataset.tag for image in images]
    inputs['night'] = [image.night for image in images]
    inputs['image_hash'] = [image_hash(image.header) for image in images]
    inputs['wcs_hash'] = [solution_hash(image.header) for image in images]
    inputs['exptime'] = [image.exptime for image in images]
    inputs['airmass'] = [-1. if image.airmass is None else image.airmass for image in images]
    fwhm_in = []
    for image in images:
        phot = image.photometry(require=False)
        fwhm_in.append(_psf(image)[1] if phot is not None and 'psf' in phot.meta else -1.)
    inputs['fwhm_in'] = fwhm_in
    inputs['fwhm_out'] = achieved
    inputs['kernel'] = kernels
    inputs['match_ee'] = ee_ratios
    inputs['offset'] = offset
    inputs['offset_err'] = offset_err
    inputs['scale'] = 10**(0.4 * offset)
    inputs['scaled'] = scaled
    inputs['nused'] = nused
    inputs['reference'] = np.arange(len(images)) == ref
    inputs['weight'] = np.array(weights) / np.max(weights)
    inputs['level'] = levels
    formatted = [
        'airmass', 'fwhm_in', 'fwhm_out', 'match_ee', 'offset', 'offset_err', 'scale', 'weight',
        'level'
    ]
    for name in formatted:
        inputs[name].info.format = '.4f'
    hdus = _product(
        images, mean, var, coverage, expmap, nreject, wcs_out, target, ref, nrej, settings,
        inputs
    )
    return hdus, inputs


def _achieved_fwhm(image, data, good, target_pix, target_beta):
    """
    Return the FWHM (arcsec) of the PSF of a matched image, or -1 if it
    cannot be measured.

    A circular Moffat model is fit to the stack of its bright, isolated stars
    (see :func:`~pipnick.photometry.psf_photometry.build_moffat_psf_model`),
    and refit with the power index fixed at that of the target, because the
    FWHM and power index of a free fit are degenerate.
    """
    try:
        fit, _ = build_moffat_psf_model(
            np.ma.MaskedArray(data, mask=np.logical_not(good)), fwhm=target_pix
        )
        # Start from the free fit, just inside the bounds of the parameters
        # (the free fit can end on a bound)
        lb, ub = fit.default_bounds()
        margin = 1e-6 * (ub - lb)
        par = fit.par.copy()
        par[4] = target_beta
        par = np.clip(par, lb + margin, ub - margin)
        fit.fit(p0=par, fix=np.array([False, False, False, False, True, False]))
    except (ValueError, RuntimeError) as e:
        log.warning(f'Could not measure the matched PSF of {image.name}: {e}')
        return -1.
    return PSFModel.from_fit(fit).fwhm * image.pixel_scale


def _check_overlap(images, reprojected):
    """
    Warn if an image does not overlap any other image of its mosaic.
    """
    if len(images) == 1:
        return
    for k, ((bk, _, _, gk), image) in enumerate(zip(reprojected, images)):
        overlaps = False
        for m, (bm, _, _, gm) in enumerate(reprojected):
            if m == k:
                continue
            rows = slice(max(bk[0].start, bm[0].start), min(bk[0].stop, bm[0].stop))
            cols = slice(max(bk[1].start, bm[1].start), min(bk[1].stop, bm[1].stop))
            if rows.stop <= rows.start or cols.stop <= cols.start:
                continue
            a = gk[
                rows.start - bk[0].start:rows.stop - bk[0].start,
                cols.start - bk[1].start:cols.stop - bk[1].start
            ]
            b = gm[
                rows.start - bm[0].start:rows.stop - bm[0].start,
                cols.start - bm[1].start:cols.stop - bm[1].start
            ]
            if np.any(a & b):
                overlaps = True
                break
        if not overlaps:
            log.warning(f'{image.name} does not overlap any other image of its mosaic.')


def _product(
    images, mean, var, coverage, expmap, nreject, wcs_out, target, ref, nrej, settings, inputs
):
    """
    Construct the FITS file of a mosaic; see :func:`mosaic_all`.
    """
    hdr = wcs_out.to_header()
    filters = sorted({str(image.filter) for image in images})
    hdr['BUNIT'] = 'electron/s'
    hdr['FILTER'] = (filters[0] if len(filters) == 1 else 'mixed', 'Filter of the images')
    hdr['NCOMBINE'] = (len(images), 'Number of images combined')
    hdr['PNVERS'] = (version, 'pipnick version (mosaic)')
    hdr['PNMOSAIC'] = (mosaic_label(images[0]), 'Mosaic label')
    hdr['PNWCSGRP'] = (wcs_group(images[0]), 'WCS group label')
    hdr['PNNIGHT'] = (images[0].night, 'UT date of the observations')
    hdr['PNPSFMAT'] = (
        settings['psf_match'] if target is not None else 'none', 'PSF matching: mosaic/group/none'
    )
    if target is not None:
        hdr['PNTGFWHM'] = (target[0], 'Target PSF FWHM (arcsec)')
        hdr['PNTGBETA'] = (target[1], 'Target PSF Moffat power index')
        hdr['PNKERNEL'] = (settings['kernel'], 'PSF matching kernel')
    hdr['PNSCALE'] = (settings['scale'], 'Flux scaling: sources, calib, or none')
    hdr['PNREFIMG'] = (images[ref].name, 'Reference image of the flux scale')
    hdr['PNNSCALE'] = (int(np.sum(inputs['scaled'])), 'Number of images on the flux scale')
    hdr['PNBKG'] = (settings['background'], 'Background subtraction')
    hdr['PNBKGMAT'] = (settings['match_background'], 'Background levels matched in overlaps')
    hdr['PNREPROJ'] = (settings['reproject'], 'Reprojection method')
    hdr['PNREJ'] = (-1. if settings['rej'] is None else settings['rej'], 'Rejection threshold')
    hdr['PNREJFRC'] = (settings['rej_frac'], 'Fractional error in the rejection')
    hdr['PNNREJ'] = (nrej, 'Number of rejected pixels')
    hdr['PNUNCTYP'] = ('density', 'UNCERT is a noise density; see pipnick docs')
    zp = images[ref].header.get('PNCALZP')
    if zp is not None:
        hdr['PHOTZP'] = (zp, 'Zero point of the reference image (e-/s), color 0')
    if len(filters) > 1:
        hdr['PNWARN'] = 'Images with different filters were combined'
    nunscaled = len(images) - int(np.sum(inputs['scaled']))
    if settings['scale'] != 'none' and nunscaled > 0:
        hdr['PNWARN2'] = f'{nunscaled} images could not be scaled'
    covered = coverage > 0
    primary = fits.PrimaryHDU(data=np.where(covered, mean, 0.).astype(np.float32), header=hdr)
    return fits.HDUList(
        [
            primary,
            fits.ImageHDU(
                data=np.sqrt(var).astype(np.float32), name='UNCERT',
                header=fits.Header({'UTYPE': 'StdDevUncertainty'})
            ),
            fits.ImageHDU(data=np.logical_not(covered).astype(np.uint8), name='MASK'),
            fits.ImageHDU(data=expmap.astype(np.float32), name='EXPMAP'),
            fits.ImageHDU(data=coverage.astype(np.int16), name='COVERAGE'),
            fits.ImageHDU(data=nreject.astype(np.int16), name='NREJECT'),
            fits.BinTableHDU(inputs, name='INPUTS')
        ]
    )


def mosaic_all(
    tables, outdir=None, mosaics=None, psf_match='mosaic', kernel='moffat', psf_margin=1.05,
    scale='sources', reference=None, match_radius=1., min_common=3, min_snr=20.,
    background='constant', match_background=False, reproject='adaptive', pixel_scale=None,
    rej=3., rej_frac=0.1, overwrite=False, qa=False, report=None
):
    r"""
    Construct the mosaics defined in the metadata tables.

    The mosaics are defined by the ``mosaic`` and ``wcs_group`` columns of the
    per-setup metadata tables; see :data:`~pipnick.metadata.mosaic_columns`.
    For each WCS group, the output WCS is found (see :func:`group_wcs`), and
    each of its mosaics is constructed (see :func:`build_mosaic`):

        1. Each image is converted to electrons per second, and its
           background is subtracted (``background``).
        2. The images are convolved to a common PSF (``psf_match``): a
           circular Moffat function with ``psf_margin`` times the largest FWHM
           of the images of the mosaic (or, for ``'group'``, of the WCS
           group), and the median power index; see
           :mod:`pipnick.mosaic.psfmatch`.  With ``'mosaic'``, the image of a
           mosaic with a single image is not convolved.
        3. The images are scaled to the flux scale of the reference image
           (``scale``): by default, from the magnitudes of the sources common
           to overlapping images, measured in apertures on the matched images;
           see :func:`measure_offsets`.
        4. The images are reprojected onto the output WCS and combined, with
           one weight per image (the inverse of its median variance) and
           outlier rejection; see :mod:`pipnick.mosaic.coadd`.

    Each mosaic is written to ``<mosaic>.fits`` in the output directory, with
    the extensions: the primary (the mosaic, in electrons per second, 0
    where there are no data), ``UNCERT`` (the error), ``MASK`` (1 where there
    are no data), ``EXPMAP`` (the total exposure time), ``COVERAGE`` (the
    number of images), ``NREJECT`` (the number of images whose pixels were
    rejected as outliers), and ``INPUTS`` (a table of the input images, with
    their provenance, PSF FWHM before and after the matching, kernel, scale,
    and weight).  The ``UNCERT`` extension is a noise *density*: the
    convolution and reprojection correlate the noise of neighboring pixels, so
    it gives the right error for sums over apertures larger than the PSF, and
    overestimates the error of single pixels.  The output WCS of each WCS group
    is written to ``<wcs_group>_wcs.fits``.  A summary is written to
    ``summary.ecsv``.

    When rerun, a mosaic is rebuilt only if its images, the settings, or the
    membership of its WCS group have changed, or ``overwrite`` is set.

    Parameters
    ----------
    tables : list
        The per-setup metadata files (see
        :func:`~pipnick.pipelines.reduction.reduce_all`).  The images of all of
        the tables with the same labels are combined.
    outdir : str, :class:`pathlib.Path`, optional
        The output directory.  If None, the ``mosaic`` subdirectory of the
        reduction directory is used, which requires a single table.
    mosaics : list, optional
        Only construct these mosaics (labels).
    psf_match : str, optional
        Match the PSFs of the images of each ``'mosaic'``, of each WCS
        ``'group'``, or not at all (``'none'``).
    kernel : str, optional
        The matching kernel: ``'moffat'`` (the kernel that matches the Moffat
        models, falling back to a Gaussian if it fails) or ``'gaussian'``.
    psf_margin : float, optional
        The target FWHM in units of the largest FWHM.
    scale : str, optional
        The flux scaling: ``'sources'``, ``'calib'`` (from the photometric
        calibration, ``PNCALZP``), or ``'none'``.
    reference : str, optional
        The reference image of the flux scale (file name or stem) of the
        mosaic that includes it.  The other mosaics, or all mosaics if None,
        use the image with the lowest airmass.
    match_radius : float, optional
        The radius (arcsec) for matching sources between images.
    min_common : int, optional
        The minimum number of sources common to a pair of images for it to
        be used for the scaling.
    min_snr : float, optional
        The minimum signal-to-noise ratio of the sources used for the scaling.
    background : str, optional
        The background subtraction: ``'constant'``, ``'2d'``, or ``'none'``.
    match_background : bool, optional
        After the scaling, match the background levels of the images: the
        additive offsets are solved from the median differences in their
        overlaps (see :func:`~pipnick.mosaic.coadd.overlap_differences`),
        with the reference image unchanged.
    reproject : str, optional
        The reprojection method: ``'adaptive'``, ``'interp'``, or ``'exact'``.
    pixel_scale : float, optional
        The output pixel scale (arcsec).  If None, the finest of the images.
    rej : float, optional
        The outlier rejection threshold; None for no rejection.
    rej_frac : float, optional
        The fractional error added to the rejection threshold, so that small
        differences between the images do not reject the cores of bright
        sources; see :func:`~pipnick.mosaic.coadd.combine`.
    overwrite : bool, optional
        Rebuild the mosaics even if they are up to date.
    qa : bool, optional
        Write a QA figure for each mosaic to the ``qa`` subdirectory.
    report : :class:`~pipnick.pkg.report.StepReport`, optional
        Collects the failures of the mosaics for the closing verdict.  A
        mosaic that fails is recorded, and the others are still constructed.
        If None, a new one is used.

    Returns
    -------
    list
        The paths to the mosaics.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if an option is invalid, or several tables are given without an
        output directory.
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the mosaics are invalid; see :func:`collect_images`.
    """
    _tables = [tables] if isinstance(tables, (str, Path)) else list(tables)
    for name, value, options in [
        ('psf_match', psf_match, ['mosaic', 'group', 'none']),
        ('kernel', kernel, ['moffat', 'gaussian']),
        ('scale', scale, ['sources', 'calib', 'none']),
        ('background', background, ['constant', '2d', 'none']),
        ('reproject', reproject, ['adaptive', 'interp', 'exact'])
    ]:
        if value not in options:
            raise PipnickParameterError(f'Invalid {name} {value}; options are {options}.')
    if outdir is None and len(_tables) > 1:
        raise PipnickParameterError('An output directory is required to combine several tables.')
    settings = {
        'psf_match': psf_match, 'kernel': kernel, 'psf_margin': float(psf_margin),
        'scale': scale, 'reference': reference, 'match_radius': float(match_radius),
        'min_common': int(min_common), 'min_snr': float(min_snr), 'background': background,
        'match_background': bool(match_background),
        'reproject': reproject, 'pixel_scale': pixel_scale,
        'rej': None if rej is None else float(rej), 'rej_frac': float(rej_frac)
    }
    _report = StepReport('MOSAICS') if report is None else report
    images = collect_images(
        _tables, require_photometry=psf_match != 'none' or scale == 'sources', report=_report
    )
    if len(images) == 0:
        _report.attention('No mosaics are defined; see pipnick_metadata --assign_mosaics.')
        _report.finish('no mosaics.')
        return []
    if outdir is None:
        _outdir = step_dir(images[0].dataset.rdxdir, 'mosaic')
    else:
        _outdir = Path(outdir).absolute()
    if not _outdir.parent.is_dir():
        raise PipnickPathError(f'{_outdir.parent} does not exist.')
    mosaic_images = _validate(images)
    if reference is not None and not any(
        reference in [image.name, image.path.stem] for image in images
    ):
        raise PipnickParameterError(f'The reference image {reference} is not in any mosaic.')
    if mosaics is not None:
        unknown = sorted(set(mosaics) - set(mosaic_images))
        if len(unknown) > 0:
            raise PipnickParameterError(f'Unknown mosaics: {", ".join(unknown)}.')

    outputs = []
    summary = []
    groups = sorted({wcs_group(image) for image in images})
    ntotal = len([m for m in mosaic_images if mosaics is None or m in mosaics])
    k = 0
    for group in groups:
        members = [image for image in images if wcs_group(image) == group]
        labels = sorted({mosaic_label(image) for image in members})
        if mosaics is not None:
            labels = [label for label in labels if label in mosaics]
        if len(labels) == 0:
            continue
        try:
            wcs_out, shape_out, ghash = group_wcs(
                members, _outdir / f'{group}_wcs.fits', pixel_scale=pixel_scale
            )
            group_target = (
                psf_target(members, margin=psf_margin) if psf_match == 'group' else None
            )
        except Exception as e:
            log.debug(f'Traceback of the failure of WCS group {group}:', exc_info=True)
            for label in labels:
                k += 1
                _report.fail(
                    f'[{k}/{ntotal}] Mosaic {label} failed: the WCS group {group} failed: '
                    f'{type(e).__name__}: {e}'
                )
                summary.append(
                    {
                        'mosaic': label, 'wcs_group': group, 'filter': '-', 'nimages': -1,
                        'nscaled': -1, 'target_fwhm': -1., 'status': 'failed',
                        'file': f'{label}.fits', 'reason': f'{type(e).__name__}: {e}'
                    }
                )
            continue
        for label in labels:
            k += 1
            progress = f'[{k}/{ntotal}]'
            ofile = _outdir / f'{label}.fits'
            these = [image for image in members if mosaic_label(image) == label]
            row = {
                'mosaic': label, 'wcs_group': group, 'filter': str(these[0].filter),
                'nimages': len(these), 'nscaled': -1, 'target_fwhm': -1., 'status': 'failed',
                'file': ofile.name, 'reason': '-'
            }
            summary.append(row)
            try:
                if psf_match == 'group':
                    target = group_target
                elif psf_match == 'mosaic' and len(these) > 1:
                    target = psf_target(these, margin=psf_margin)
                else:
                    # A single image is not convolved to match itself; it is
                    # only reprojected
                    target = None
                mhash = json_hash(
                    {
                        'settings': settings, 'group': ghash, 'target': target,
                        'images': sorted(_image_key(image) for image in these)
                    }
                )
                if ofile.is_file() and not overwrite:
                    hdr = fits.getheader(ofile)
                    if hdr.get('PNMOSHSH') == mhash:
                        log.info(f'{progress} {ofile.name} is up to date.')
                        _report.check_version(ofile.name, hdr.get('PNVERS'))
                        row.update(_summary_values(ofile))
                        row['status'] = 'up to date'
                        outputs.append(ofile)
                        continue
                log.info(f'{progress} Constructing mosaic {label} from {len(these)} images.')
                hdus, inputs = build_mosaic(these, wcs_out, shape_out, target, settings)
            except Exception as e:
                log.debug(f'Traceback of the failure of mosaic {label}:', exc_info=True)
                _report.fail(
                    f'Mosaic {label} failed: {type(e).__name__}: {e}.  Continuing with the other '
                    'mosaics.'
                )
                row['reason'] = f'{type(e).__name__}: {e}'
                continue
            hdus[0].header['PNMOSHSH'] = (mhash, 'Hash of the mosaic inputs and settings')
            hdus[0].header['PNGRPHSH'] = (ghash, 'Hash of the WCS group membership')
            _outdir.mkdir(parents=True, exist_ok=True)
            hdus.writeto(ofile, overwrite=True)
            log.info(f'Wrote {ofile.parent.name}/{ofile.name}')
            row.update(_summary_values(ofile))
            row['status'] = 'constructed'
            outputs.append(ofile)
            if qa:
                _write_qa(hdus, inputs, _outdir / 'qa' / f'{label}.png')
    _summarize(summary, _outdir)
    nbuilt = sum(row['status'] == 'constructed' for row in summary)
    nup = sum(row['status'] == 'up to date' for row in summary)
    nfail = sum(row['status'] == 'failed' for row in summary)
    _report.finish(
        f'{len(summary)} mosaics ({nbuilt} constructed, {nup} up to date, {nfail} failed); '
        f'summary in {_outdir.name}/summary.ecsv.',
        next_step='pipnick_rgb, to combine the mosaics of a WCS group into a colour image'
    )
    return outputs


def _image_key(image):
    """
    Return the provenance of an image used in the hash of its mosaic: its
    path, the hashes of its content and astrometric solution, and the hash of
    its photometry settings.
    """
    phot = image.photometry(require=False)
    return [
        str(image.path), image_hash(image.header), solution_hash(image.header),
        None if phot is None else phot.meta.get('settings_hash')
    ]


def _summary_values(ofile):
    """
    Return the summary values of a mosaic file; see :func:`_summarize`.
    """
    hdr = fits.getheader(ofile)
    return {
        'nscaled': int(hdr.get('PNNSCALE', -1)), 'target_fwhm': float(hdr.get('PNTGFWHM', -1.))
    }


def _summarize(summary, outdir):
    """
    Log the summary of the mosaics, and write it to ``summary.ecsv``.
    """
    if len(summary) == 0:
        return
    tbl = Table(rows=summary)
    for name in ['nscaled', 'target_fwhm']:
        tbl[name] = np.ma.MaskedArray(tbl[name], mask=np.asarray(tbl[name]) < 0)
    tbl['target_fwhm'].info.format = '.2f'
    tbl.meta['pipnick_version'] = version
    outdir.mkdir(parents=True, exist_ok=True)
    tbl.write(outdir / 'summary.ecsv', format='ascii.ecsv', overwrite=True)
    log.info(
        f'Mosaic summary (also in {outdir / "summary.ecsv"}):\n'
        + '\n'.join(tbl.pformat(max_lines=-1, max_width=-1))
    )


def _write_qa(hdus, inputs, ofile):
    """
    Write the QA figure of a mosaic: the mosaic, its coverage, the magnitude
    offsets of the images, and their PSF FWHM before and after the matching.
    """
    data = hdus[0].data
    good = hdus['MASK'].data == 0
    fig = Figure(figsize=(12, 10))
    axes = fig.subplots(2, 2)
    vmin, vmax = ZScaleInterval().get_limits(data[good]) if np.any(good) else (0., 1.)
    axes[0,0].imshow(data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
    axes[0,0].set_title(f'{hdus[0].header["PNMOSAIC"]} ({hdus[0].header["FILTER"]})')
    im = axes[0,1].imshow(hdus['COVERAGE'].data, origin='lower', cmap='viridis')
    fig.colorbar(im, ax=axes[0,1], label='number of images')
    axes[0,1].set_title('coverage')
    index = np.arange(len(inputs))
    axes[1,0].errorbar(index, inputs['offset'], yerr=inputs['offset_err'], fmt='o')
    unscaled = np.logical_not(np.asarray(inputs['scaled']))
    axes[1,0].scatter(
        index[unscaled], np.asarray(inputs['offset'])[unscaled], marker='x', color='C3', s=60,
        label='not scaled'
    )
    axes[1,0].set_xlabel('image')
    axes[1,0].set_ylabel('magnitude offset (mag)')
    axes[1,0].set_title(f'reference: {hdus[0].header["PNREFIMG"]}', fontsize=9)
    axes[1,1].plot(index, inputs['fwhm_in'], 'o', label='input')
    fwhm_out = np.asarray(inputs['fwhm_out'])
    axes[1,1].plot(index[fwhm_out > 0], fwhm_out[fwhm_out > 0], 's', label='matched')
    if 'PNTGFWHM' in hdus[0].header:
        axes[1,1].axhline(hdus[0].header['PNTGFWHM'], color='k', lw=0.5, label='target')
    axes[1,1].set_xlabel('image')
    axes[1,1].set_ylabel('PSF FWHM (arcsec)')
    axes[1,1].legend()
    fig.tight_layout()
    ofile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ofile)
    log.info(f'Wrote {ofile}')
