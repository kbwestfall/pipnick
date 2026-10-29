"""
Point-spread-function (PSF) photometry.

The PSF of an image is modeled as a 2D Moffat function fit to a stack of
bright, isolated sources (:func:`build_moffat_psf_model`), and the sources in
the image are measured by fitting this PSF with ``photutils``
(:func:`psf_analysis`).
"""

import numpy as np

from astropy.modeling.fitting import TRFLSQFitter
from photutils.background import LocalBackground, MADStdBackgroundRMS, MMMBackground
from photutils.psf import IterativePSFPhotometry, SourceGrouper

from pipnick import log
from pipnick.datamodel.psf import PSFModel
from pipnick.image.background import subtract_background
from pipnick.image.detect import SourceFinder, find_sources
from pipnick.photometry import starfind
from pipnick.photometry.moffat import FitEllipticalMoffat2D, FitMoffat2D


def build_moffat_psf_model(
    image, thresh=15.0, fwhm=7., fittype='circ', stack='mean', pdf_file=None
):
    """
    Use well-detected and well-separated sources to model the image PSF as a
    2D Moffat function.

    Sources are detected with :func:`~pipnick.image.detect.find_sources`
    (with the global background, and the masked pixels set to the background
    level) and selected with
    :func:`~pipnick.photometry.starfind.select_psf_stars`, background-subtracted
    stamps are extracted around them
    (:func:`~pipnick.photometry.starfind.extract_stamps`), and the stamps with
    less than 1% masked pixels are stacked and fit.

    Parameters
    ----------
    image : str, :class:`~pathlib.Path`, :class:`~astropy.nddata.CCDData`, array-like
        Image for PSF analysis; see
        :func:`~pipnick.photometry.starfind.as_masked_array`.
    thresh : float, optional
        Detection threshold in units of background standard deviation.
    fwhm : float, optional
        Initial estimate of the FWHM in pixels.
    fittype : str, optional
        Type of PSF fit to use: ``'circ'`` for circular or ``'ellip'`` for
        elliptical.
    stack : str, optional
        How the stamps are stacked: ``'mean'`` or ``'median'``.
    pdf_file : str, :class:`pathlib.Path`, optional
        If provided, plot the stacked stamps to this file; see
        :func:`~pipnick.photometry.starfind.plot_stamps`.

    Returns
    -------
    fit : object
        The fit to the stack, a :class:`~pipnick.photometry.moffat.FitMoffat2D`
        or :class:`~pipnick.photometry.moffat.FitEllipticalMoffat2D` object.
        The stack, best-fit parameters, and model are ``fit.c``, ``fit.par``,
        and ``fit.model()``.
    nstack : int
        The number of stacked sources.

    Raises
    ------
    ValueError
        Raised if ``fittype`` or ``stack`` is not valid, or if no sources are
        good enough to model the PSF.
    """
    if fittype not in ['circ', 'ellip']:
        raise ValueError(f'Unknown fit type {fittype}; must be circ or ellip.')
    if stack not in ['mean', 'median']:
        raise ValueError(f'Unknown stacking method {stack}; must be mean or median.')
    _image = starfind.as_masked_array(image)
    sources = find_sources(_image, fwhm=fwhm, thresh=thresh)
    source_data = starfind.select_psf_stars(_image, sources, fwhm=fwhm)
    log.debug(f'Found {len(source_data)} well-separated sources with well-detected fluxes.')

    log.debug('Extracting image stamps.')
    _, stamps = starfind.extract_stamps(_image, source_data)
    npix_stamp = np.prod(stamps.shape[1:])
    mask_fraction = np.sum(np.ma.getmaskarray(stamps), axis=(1,2)) / npix_stamp
    indx = source_data['good'].data & (mask_fraction < 0.01)
    nstack = int(np.sum(indx))
    log.debug(f'Found {nstack} good sources to stack for the PSF measurement.')
    if nstack == 0:
        raise ValueError('No sources are good enough to model the PSF.')

    if pdf_file is not None:
        starfind.plot_stamps(stamps[indx], pdf_file, vmin=None, vmax=None)

    if stack == 'mean':
        stacked = np.ma.mean(stamps[indx], axis=0)
    else:
        stacked = np.ma.median(stamps[indx], axis=0)
    fit = FitMoffat2D(stacked) if fittype == 'circ' else FitEllipticalMoffat2D(stacked)
    fit.fit()
    return fit, nstack


def psf_analysis(
    image, error=None, thresh=10.0, guess_fwhm=7., mode='all', fittype='circ', stack='mean',
    min_sep=1., group_sep=2., psf_thresh=15.
):
    """
    Perform PSF photometry on an image, including source detection and PSF
    fitting.

    The PSF is modeled with :func:`build_moffat_psf_model`.  A smooth
    background is subtracted (see
    :func:`~pipnick.image.background.subtract_background`), and the sources
    are found (with :class:`~pipnick.image.detect.SourceFinder`, above
    ``thresh`` times a robust estimate of the noise of the subtracted image)
    and fit with :class:`photutils.psf.IterativePSFPhotometry`, with the PSF
    shape fixed and the flux and position of each source free.  The residual
    images of the later iterations are searched with the same threshold.
    Sources closer than ``group_sep`` FWHM are fit simultaneously.

    Parameters
    ----------
    image : str, :class:`~pathlib.Path`, :class:`~astropy.nddata.CCDData`, array-like
        Image for PSF analysis; see
        :func:`~pipnick.photometry.starfind.as_masked_array`.
    error : :class:`numpy.ndarray`, optional
        The 1-sigma error in each pixel, used to weight the fits and to
        compute the flux errors.
    thresh : float, optional
        Detection threshold in units of the background standard deviation.
    guess_fwhm : float, optional
        Initial guess for the FWHM of the PSF in pixels.
    mode : str, optional
        Mode for the iterative PSF fitting (``'all'`` or ``'new'``); see
        :class:`photutils.psf.IterativePSFPhotometry`.
    fittype : str, optional
        Type of PSF fit to use: ``'circ'`` for circular or ``'ellip'`` for
        elliptical.
    stack : str, optional
        How the stamps used to model the PSF are stacked; see
        :func:`build_moffat_psf_model`.
    min_sep : float, optional
        The minimum separation of detected sources, in units of the FWHM.
    group_sep : float, optional
        Sources closer than this, in units of the FWHM, are fit
        simultaneously; see :class:`photutils.psf.SourceGrouper`.
    psf_thresh : float, optional
        Detection threshold, in units of the background standard deviation,
        for the sources used to model the PSF.

    Returns
    -------
    fit : object
        The fit to the PSF, a :class:`~pipnick.photometry.moffat.FitMoffat2D`
        or :class:`~pipnick.photometry.moffat.FitEllipticalMoffat2D` object;
        see :func:`build_moffat_psf_model`.
    nstack : int
        The number of sources used to model the PSF.
    fwhm : float
        The FWHM of the PSF in pixels; see
        :attr:`~pipnick.datamodel.psf.PSFModel.fwhm`.
    phot_data : :class:`astropy.table.Table`
        The results of the PSF photometry; see
        :class:`photutils.psf.IterativePSFPhotometry`.  Failed fits have NaN
        values, which callers must mask.
    data : :class:`numpy.ndarray`
        The background-subtracted image used for the photometry, with masked
        pixels set to 0.
    """
    _image = starfind.as_masked_array(image)
    mask = np.ma.getmaskarray(_image)
    fit, nstack = build_moffat_psf_model(
        _image, thresh=psf_thresh, fwhm=guess_fwhm, fittype=fittype, stack=stack
    )
    psf = PSFModel.from_fit(fit)
    fwhm = psf.fwhm
    log.debug(f'Measured Moffat PSF FWHM = {fwhm:.2f} pixels from {nstack} sources.')

    # Scales set by the FWHM
    aper_size = fwhm * 1.8
    local_bkg_range = (3 * fwhm, 6 * fwhm)
    win = int(np.ceil(2 * fwhm))
    if win % 2 == 0:
        win += 1  # Ensure window size is odd

    data = subtract_background(_image.filled(0.), mask=mask)
    # A robust estimate of the background noise
    std = MADStdBackgroundRMS()(np.ma.MaskedArray(data, mask=mask))

    # NOTE: DAOStarFinder (see find_sources) applies a matched filter before
    # thresholding; IRAFStarFinder does not, and so detects many more
    # spurious sources
    finder = SourceFinder(fwhm, thresh * std, min_sep=min_sep * fwhm)
    grouper = SourceGrouper(min_separation=group_sep * fwhm)
    local_bkg = LocalBackground(*local_bkg_range, bkg_estimator=MMMBackground())

    # The PSF model, with its shape fixed by the fit to the stack
    phot = IterativePSFPhotometry(
        psf.to_photutils(), win, finder, grouper=grouper, local_bkg_estimator=local_bkg,
        fitter=TRFLSQFitter(), aperture_radius=aper_size, mode=mode
    )
    phot_data = phot(data, mask=mask, error=error)
    return fit, nstack, fwhm, phot_data, data


def consolidate_groups(phot_data, preserve=None):
    """
    Consolidate groups of sources in photometric source into one source.

    This function combines data for sources that are identified as belonging
    to the same group by summing fluxes and weighted-averaging other parameters.

    Parameters
    ----------
    phot_data : :class:`astropy.table.Table`
        Table containing photometric data of detected sources; see
        :func:`psf_analysis`.
    preserve : list, optional
        List of group IDs to preserve during consolidation.

    Returns
    -------
    :class:`astropy.table.Table`
        A new table containing consolidated photometric data.
    """
    _preserve = [] if preserve is None else preserve
    group_data = phot_data[phot_data['group_size'] > 1]
    # Exclude preserved groups
    group_data = group_data[np.logical_not(np.isin(group_data['group_id'], _preserve))]

    group_ids = sorted(set(group_data['group_id']))
    log.debug(f'Consolidating {len(group_ids)} groups of blended sources: {group_ids}')

    # Start with non-group data
    new_data = phot_data[np.logical_not(np.isin(phot_data['group_id'], group_ids))]

    # Consolidate data for each group
    weighted = [
        'local_bkg', 'x_init', 'y_init', 'flux_init', 'x_fit', 'y_fit', 'qfit', 'cfit',
        'reduced_chi2'
    ]
    summed = ['flux_fit', 'x_err', 'y_err', 'flux_err', 'n_pixels_fit']
    for id in group_ids:
        group = group_data[group_data['group_id'] == id]
        weights = np.absolute(group['flux_fit'])
        new_row = {'id': id, 'group_id': id, 'group_size': group['group_size'][0]}
        for name in group.colnames:
            if name in weighted:
                new_row[name] = np.average(group[name], weights=weights)
            elif name in summed:
                new_row[name] = np.sum(group[name])
            elif name == 'iter_detected':
                new_row[name] = group[name][-1]
            elif name == 'flags':
                new_row[name] = np.bitwise_or.reduce(group[name])
        new_data.add_row(new_row)

    new_data.sort('group_id')
    return new_data
