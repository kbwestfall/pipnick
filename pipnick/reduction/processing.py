"""
Process raw frames and combine them.

Each frame is overscan-subtracted, trimmed, bias-subtracted, field-flattened,
and cleaned of cosmic rays, with its uncertainty propagated (see
:func:`process_frame`), and frames are combined with a sigma-clipped mean (see
:func:`combine_frames`).  The camera properties (sections, bad pixels, gain,
read noise, and saturation) come from the camera configuration; see
:mod:`pipnick.cameras`.
"""

from pathlib import Path
import warnings

import numpy as np

from astropy.nddata import CCDData, StdDevUncertainty
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy import units
import ccdproc

from pipnick import cameras
from pipnick import log
from pipnick.pkg.exceptions import PipnickDataModelError
from pipnick.reduction.combiner import Combiner


def init_ccddata(frame, bpm=None, saturation=None, gain=None, readnoise=None):
    """
    Read data from a provided file.
    
    This function also uses the header information to identify the camera used
    to acquire the data.  If the camera is successfully identified, the
    camera information is used to 

        - apply any known bad-pixel mask
        - mask saturated pixels
        - convert the units from DN to counts by multiplying by the camera gain

    The relevant keyword arguments can be used to supplement information not
    available for a given camera or when the camera cannot be identified.  Note
    that values provided via a keyword *always* take precedence.

    Parameters
    ----------
    frame : str, :class:`pathlib.Path`
        Path to the FITS file.

    Returns
    -------
    :class:`astropy.nddata.CCDData`
        Initialized and processed CCDData object.
    """
    _frame = Path(frame).absolute()

    log.debug(f'Initializing CCD imaging data in {_frame.name}.')
    ccd = CCDData.read(_frame, unit=units.adu)

    # Determine the camera used
    try:
        camera = cameras.identify_camera(ccd.header)
    except PipnickDataModelError:
        warnings.warn(f'Unable to determine camera used for {_frame.name}.')
        camera = None
        readout = None
    else:
        ccd.meta['camera'] = camera.name
        readout = camera.readout(ccd.header)

    # Mask bad pixels
    _bpm = camera.bpm(ccd.header) if bpm is None and camera is not None else bpm
    if _bpm is not None:
        ccd.mask = _bpm

    # Mask saturated pixels; saturation is in ADU
    _sat = readout.saturation if saturation is None and readout is not None else saturation
    if _sat is not None:
        if ccd.mask is None:
            ccd.mask = np.zeros(ccd.data.shape, dtype=bool)
        saturated = ccd.data >= _sat
        nsat = int(np.sum(saturated))
        # The number of saturated pixels that are not known bad pixels, which
        # can read as saturated; see process_frame
        ccd.meta['PNNSAT'] = int(np.sum(saturated & np.logical_not(ccd.mask)))
        if nsat > 0:
            log.debug(
                f'{_frame.name} has {nsat} saturated pixels (>= {_sat} ADU), '
                f'{nsat - ccd.meta["PNNSAT"]} of which are known bad pixels.'
            )
        ccd.mask[saturated] = True

    # Get and save the RN (in electrons) and gain (in electrons/ADU)
    ccd.meta['PROCGAIN'] = readout.gain if gain is None and readout is not None else gain
    ccd.meta['PROCRN'] = (
        readout.readnoise if readnoise is None and readout is not None else readnoise
    )

    # Apply gain if available
    if ccd.meta['PROCGAIN'] is not None:
        ccd.data = ccd.data * ccd.meta['PROCGAIN']
        ccd.unit *= (units.electron / units.adu)

    # Done
    return ccd



def process_frame(raw_frame, bias=None, flat=None, flag_cosmics=True, photon_noise=True):
    r"""
    Perform basic image processing

    This includes overscan subtraction, trimming, bias subtraction, and
    field-flattening.  One line is logged for the frame, with the mean
    overscan level (in ADU), the number of saturated pixels (that are not
    known bad pixels), the calibrations applied, and the number of pixels
    affected by cosmic rays.  The overscan level is saved to the ``PNOSCAN``
    header keyword, and the numbers of saturated and cosmic-ray pixels to the
    ``PNNSAT`` and ``PNNCR`` keywords of the frame's header (which are not
    written to the combined output; see :func:`~pipnick.pipelines.reduction.reduce_all`).

    The uncertainty of each pixel is propagated through the processing:

        - After the overscan subtraction, the variance is the read noise,
          :math:`\sigma_{\rm RN}^2`, plus the error in the subtracted overscan
          level, :math:`\sigma_{\rm RN}^2 / N_{\rm os}`.  The overscan level
          of each row is the mean of its :math:`N_{\rm os}` overscan pixels
          (see :func:`ccdproc.subtract_overscan`), so this error is the same
          for all the pixels in a row.  Without an overscan region, the
          variance is :math:`\sigma_{\rm RN}^2`.
        - The bias subtraction adds the variance of the bias image.
        - If ``photon_noise`` is True, the Poisson variance of the
          bias-subtracted counts (in electrons, and not less than 0) is
          added.
        - The field-flattening propagates the uncertainties of the image and
          the flat.

    The uncertainty is not set if the read noise is unknown.

    Parameters
    ----------
    raw_frame : str, :class:`pathlib.Path`
        Data file to process.
    flag_cosmics : bool, optional
        Identify and mask cosmic rays
    bias : :class:`astropy.nddata.CCDData`, optional
        Image to use for the bias subtraction
    flat : :class:`astropy.nddata.CCDData`, optional
        Image to use for the field-flattening
    photon_noise : bool, optional
        Include the Poisson noise of the counts in the uncertainty.  This
        should be False for biases, which have no photons.

    Returns
    -------
    :class:`astropy.nddata.CCDData`
        Processed CCDData object with overscan subtracted and image trimmed.
    """
    ccd = init_ccddata(raw_frame)
    name = Path(raw_frame).name
    # The summary of the processing, logged at the end
    summary = []
    if ccd.meta.get('PNNSAT', 0) > 0:
        summary.append(f'{ccd.meta["PNNSAT"]} saturated pixels')
    readnoise = ccd.meta.get('PROCRN')
    # The number of overscan pixels averaged to set the overscan level of each
    # row; 0 if there is no overscan subtraction
    noverscan = 0
    if 'camera' not in ccd.meta:
        log.warning('Unknown camera; cannot determine the overscan and data sections!')
    else:
        datasec, biassec = cameras.get_camera(ccd.meta['camera']).sections(ccd.header)
        if biassec is None:
            log.warning('Frame has no overscan region.')
        else:
            overscan = ccdproc.trim_image(ccd, fits_section=biassec)
            noverscan = overscan.shape[1]
            level = _masked_mean(overscan)
            if ccd.meta['PROCGAIN'] is not None:
                level /= ccd.meta['PROCGAIN']
            ccd.meta['PNOSCAN'] = (level, 'Mean overscan level (ADU)')
            log.debug(f'Subtracting overscan; mean level is {level:.1f} ADU.')
            summary.insert(0, f'overscan {level:.1f} ADU')
            ccd = ccdproc.subtract_overscan(ccd, fits_section=biassec, overscan_axis=1)
        ccd = ccdproc.trim_image(ccd, fits_section=datasec)
        log.debug(f'Trimmed image to {ccd.shape[1]} x {ccd.shape[0]} (columns x rows).')
    if readnoise is None:
        log.warning('Unknown read noise; the uncertainty of the image is not set.')
    else:
        # Read noise, including the error in the subtracted overscan level
        variance = readnoise**2 * (1 + (0. if noverscan == 0 else 1 / noverscan))
        ccd.uncertainty = StdDevUncertainty(np.full(ccd.shape, np.sqrt(variance)))
    if bias is not None:
        log.debug('Subtracting super bias')
        summary.append('bias subtracted')
        ccd = ccdproc.subtract_bias(ccd, bias)
    if photon_noise and ccd.uncertainty is not None:
        sigma = np.asarray(ccd.uncertainty.array)
        ccd.uncertainty = StdDevUncertainty(
            np.sqrt(sigma**2 + np.clip(np.asarray(ccd.data), 0., None))
        )
    if flat is not None:
        log.debug('Applying flat-field correction')
        summary.append('flat-fielded')
        ccd = ccdproc.flat_correct(ccd, flat)
    if flag_cosmics:
        nmasked = 0 if ccd.mask is None else int(np.sum(ccd.mask))
        ccd = ccdproc.cosmicray_lacosmic(
            ccd, gain_apply=False, gain=1.0, readnoise=ccd.meta['PROCRN'], verbose=False
        )
        ncosmic = 0 if ccd.mask is None else int(np.sum(ccd.mask)) - nmasked
        ccd.meta['PNNCR'] = ncosmic
        summary.append(f'{ncosmic} cosmic-ray pixels ({100 * ncosmic / ccd.data.size:.2f}%)')
    log.info(f'{name}: ' + ('processed' if len(summary) == 0 else '; '.join(summary)) + '.')
    return ccd



def _masked_mean(ccd):
    """
    Return the mean of the unmasked pixels in an image.
    """
    data = np.ma.MaskedArray(ccd.data, mask=ccd.mask)
    return float(np.ma.mean(data))



def process_frames(raw_frames, bias=None, flat=None, flag_cosmics=True, photon_noise=True):
    """
    Perform basic image processing of a list of frames; see
    :func:`process_frame`.

    Parameters
    ----------
    raw_frames : list
        List of file names to process.
    bias : :class:`astropy.nddata.CCDData`, optional
        Image to use for the bias subtraction
    flat : :class:`astropy.nddata.CCDData`, optional
        Image to use for the field-flattening
    flag_cosmics : bool, optional
        Identify and mask cosmic rays

    Returns
    -------
    list
        The processed :class:`astropy.nddata.CCDData` objects.
    """
    # TODO: For *many* images, this may be too memory intensive.
    return [
        process_frame(
            frame, bias=bias, flat=flat, flag_cosmics=flag_cosmics, photon_noise=photon_noise
        )
        for frame in raw_frames
    ]



def combine_frames(frames, scale=None):
    """
    Combine a list of processed frames using a sigma-clipped mean.

    Outliers are iteratively rejected at 5 sigma.  The header and WCS of the
    combined frame are those of the first frame, with ``NCOMBINE`` set to the
    number of frames.  A single frame is returned without combination (but
    with any scaling applied).

    The uncertainty of the combined frame is propagated from the
    uncertainties of the (scaled) frames that are not rejected at each pixel;
    see :class:`~pipnick.reduction.combiner.Combiner`.  If any frame has no
    uncertainty, the ``ccdproc`` default is used instead: the standard
    deviation of the combined values divided by the square root of their
    number.

    Parameters
    ----------
    frames : list
        The :class:`astropy.nddata.CCDData` objects to combine.
    scale : array-like, optional
        Factors used to multiply each frame before combining them.  If None,
        the frames are not scaled.

    Returns
    -------
    :class:`astropy.nddata.CCDData`
        Combined CCDData object.
    """
    # NOTE: The arithmetic and the combination do not propagate the header
    meta = frames[0].meta.copy()
    wcs = frames[0].wcs
    _frames = frames if scale is None else [f.multiply(s) for s, f in zip(scale, frames)]

    if len(_frames) == 1:
        combined = _frames[0]
    else:
        log.debug('Identifying outliers using sigma clipping')
        combiner = Combiner(_frames)

        old_n_masked = 0
        new_n_masked = 1
        while new_n_masked > old_n_masked:
            # TODO:
            #   - does ccdproc provide an iterative sigma clipper?
            #   - allow user to change the clipping parameters
            combiner.data_arr = sigma_clip(
                combiner.data_arr, sigma=5., maxiters=1, cenfunc=np.ma.mean, stdfunc=np.ma.std,
                axis=0, masked=True, copy=False
            )
            old_n_masked = new_n_masked
            new_n_masked = combiner.data_arr.mask.sum()

        log.debug('Creating sigma-clipped average image.')
        propagate = combiner.uncertainty_arr is not None
        if not propagate:
            log.warning(
                'Not all frames have uncertainties; the uncertainty of the combined frame is '
                'the scatter of the combined values.'
            )
        combined = combiner.average_combine(
            scale_func=np.ma.mean, propagate_uncertainty=propagate
        )

    # TODO: Fill the masked regions?
    combined.meta = meta
    combined.meta['NCOMBINE'] = len(frames)
    combined.wcs = wcs
    return combined



def report_bias_stats(frames, files):
    r"""
    Report the bias level of the bias frames and the read noise measured from
    them.

    The bias level of each frame is the mean of the unmasked pixels in the
    data section before overscan subtraction, in ADU, for comparison with the
    overscan level.  The read noise is measured from the difference of each
    consecutive pair of frames as :math:`\sigma_{\rm diff}/\sqrt{2}`, where
    :math:`\sigma_{\rm diff}` is the standard deviation of the difference
    after rejecting 5-sigma outliers (e.g., cosmic rays).  The frames include
    the nominal gain, so the read noise is in electrons.

    Parameters
    ----------
    frames : list
        The processed bias frames; see :func:`process_frames`.
    files : list
        The raw files.

    Returns
    -------
    level : float
        The median bias level of the frames (ADU).
    readnoise : float
        The median read noise (e-), or None if there are fewer than two frames
        or the gain is unknown.
    """
    levels = []
    for f, frame in zip(files, frames):
        gain = frame.meta['PROCGAIN']
        residual = _masked_mean(frame) / (1. if gain is None else gain)
        oscan = frame.meta.get('PNOSCAN', 0.)
        levels.append(oscan + residual)
        log.debug(
            f'{Path(f).name}: mean bias level is {oscan + residual:.1f} ADU ({residual:.2f} ADU '
            'after overscan subtraction).'
        )
    log.info(
        f'Bias level of the {len(levels)} frames: median {np.median(levels):.1f} ADU (range '
        f'{np.min(levels):.1f}-{np.max(levels):.1f}).'
    )
    if len(frames) < 2:
        log.info('At least two bias frames are needed to measure the read noise.')
        return float(np.median(levels)), None
    if frames[0].meta['PROCGAIN'] is None:
        log.info('The gain is unknown, so the read noise cannot be measured in electrons.')
        return float(np.median(levels)), None

    readnoise = []
    for i in range(len(frames) - 1):
        mask = np.zeros(frames[i].shape, dtype=bool)
        for frame in frames[i:i+2]:
            if frame.mask is not None:
                mask |= frame.mask
        diff = np.ma.MaskedArray(frames[i+1].data - frames[i].data, mask=mask)
        _, _, std = sigma_clipped_stats(diff, sigma=5.)
        readnoise.append(std / np.sqrt(2.))
    nominal = frames[0].meta['PROCRN']
    log.info(
        f'Read noise from {len(readnoise)} pairs of bias frames: median {np.median(readnoise):.2f}'
        f' e- (range {np.min(readnoise):.2f}-{np.max(readnoise):.2f}); nominal value is '
        f'{nominal} e-.'
    )
    return float(np.median(levels)), float(np.median(readnoise))

