"""
Photometry of the reduced science images.

For each image, the point-spread function (PSF) is modeled, and the sources
are measured by PSF and aperture photometry; see
:func:`~pipnick.photometry.psf_photometry.psf_analysis`.  The results are
written to a source table for each image (``photometry/<image>_phot.ecsv`` in
the reduction directory).

The pixel coordinates of the sources are always provided.  Their celestial
coordinates are provided only if the image has an astrometric solution found
by pipnick (i.e., its header has the ``PNASTMTH`` keyword; see
:func:`~pipnick.pipelines.astrometry.astrometry_all`), and are masked
otherwise.  The WCS that the reduced images inherit from the raw headers is the
nominal telescope solution, which is not used.
"""

import numpy as np

from astropy.table import MaskedColumn, Table

from pipnick import log
from pipnick.datamodel.catalog import SourceCatalog, celestial_coordinates
from pipnick.datamodel.dataset import DataSet
from pipnick.datamodel.psf import PSFModel
from pipnick.image.io import invalid_error
from pipnick.photometry.aperture_photometry import aperture_analysis
from pipnick.photometry.psf_photometry import consolidate_groups, psf_analysis
from pipnick.pkg.exceptions import PipnickPathError
from pipnick.pkg.provenance import image_hash, json_hash, solution_hash
from pipnick.pkg.report import StepReport
from pipnick.pkg.summary import Summary
from pipnick.pkg.version import version
from pipnick.qa.psf import plot_psf
from pipnick.qa.sources import plot_sources


def _masked(values):
    """
    Return a masked array with the non-finite values (e.g., from failed fits)
    masked and replaced by 0.
    """
    _values = np.asarray(values, dtype=float)
    bad = np.logical_not(np.isfinite(_values))
    return np.ma.MaskedArray(np.where(bad, 0., _values), mask=bad)


def measure_image(image, settings):
    """
    Measure the sources in a reduced image.

    The object, filter, exposure time, and airmass recorded with the sources
    are the metadata table values of the image (see
    :attr:`~pipnick.datamodel.image.ReducedImage.values`), and the pixel scale
    is its nominal scale.

    Parameters
    ----------
    image : :class:`~pipnick.datamodel.image.ReducedImage`
        The reduced image.
    settings : dict
        The photometry settings; see :func:`photometry_all`.

    Returns
    -------
    tbl : :class:`astropy.table.Table`
        The source table; see :func:`photometry_all`.
    fit : object
        The fit to the PSF, a :class:`~pipnick.photometry.moffat.FitMoffat2D`
        or :class:`~pipnick.photometry.moffat.FitEllipticalMoffat2D` object.
    data : :class:`numpy.ndarray`
        The background-subtracted image.
    mask : :class:`numpy.ndarray`
        The image mask.
    """
    path = image.path
    values = image.values
    mjd = image.mjd
    scale = image.nominal_scale
    data, mask, error, hdr = image.read()
    # The pixel errors propagated by the reduction (see
    # pipnick.pipelines.reduction.process_frame)
    if error is None:
        log.warning(f'{path.name} has no uncertainty; the flux errors are not weighted.')
    else:
        # Pixels without a valid error are ignored
        bad = invalid_error(error)
        mask = mask | bad
        error = np.where(bad, 1., error)
    data = np.ma.MaskedArray(data, mask=mask)

    fit, nstack, fwhm, phot, sub = psf_analysis(
        data, error=error, thresh=settings['thresh'], guess_fwhm=settings['fwhm'] / scale,
        mode=settings['mode'], fittype=settings['fittype'], stack=settings['stack'],
        min_sep=settings['min_sep'], group_sep=settings['group_sep'],
        psf_thresh=settings['psf_thresh']
    )
    if settings['group']:
        phot = consolidate_groups(phot)

    phot, removed = vet_sources(phot, data.shape, fwhm)
    x = np.asarray(phot['x_fit'], dtype=float)
    y = np.asarray(phot['y_fit'], dtype=float)

    flux_aper, flux_aper_err = aperture_analysis(
        sub, x, y, settings['aper_radius'] * fwhm, local_bkg=np.asarray(phot['local_bkg']),
        error=error, mask=mask
    )
    flux_psf = _masked(phot['flux_fit'])
    flux_psf_err = _masked(phot['flux_err'])

    # Instrumental magnitudes, per second
    exptime = values['exptime']
    good = np.logical_not(np.ma.getmaskarray(flux_psf)) & (flux_psf.filled(0.) > 0)
    mag = np.ma.MaskedArray(np.zeros(len(x)), mask=np.logical_not(good))
    mag_err = np.ma.MaskedArray(np.zeros(len(x)), mask=np.logical_not(good))
    if exptime is not None and np.any(good):
        mag[good] = -2.5 * np.log10(flux_psf[good] / exptime)
        mag_err[good] = 2.5 / np.log(10.) * flux_psf_err[good] / flux_psf[good]
        mag_err[good & np.ma.getmaskarray(flux_psf_err)] = np.ma.masked

    ra, dec = celestial_coordinates(hdr, x, y)
    tbl = Table()
    tbl['id'] = np.arange(1, len(x) + 1)
    for name in ['group_id', 'group_size', 'iter_detected']:
        if name in phot.colnames:
            tbl[name] = np.asarray(phot[name])
    tbl['x'] = MaskedColumn(x, unit='pix')
    tbl['y'] = MaskedColumn(y, unit='pix')
    tbl['x_err'] = MaskedColumn(_masked(phot['x_err']), unit='pix')
    tbl['y_err'] = MaskedColumn(_masked(phot['y_err']), unit='pix')
    tbl['ra'] = ra
    tbl['dec'] = dec
    tbl['flux_psf'] = MaskedColumn(flux_psf, unit='electron')
    tbl['flux_psf_err'] = MaskedColumn(flux_psf_err, unit='electron')
    tbl['flux_aper'] = MaskedColumn(flux_aper, unit='electron')
    tbl['flux_aper_err'] = MaskedColumn(flux_aper_err, unit='electron')
    tbl['mag_inst'] = MaskedColumn(mag, unit='mag')
    tbl['mag_inst_err'] = MaskedColumn(mag_err, unit='mag')
    tbl['local_bkg'] = MaskedColumn(_masked(phot['local_bkg']), unit='electron')
    for name in ['qfit', 'cfit', 'reduced_chi2']:
        tbl[name] = MaskedColumn(_masked(phot[name]))
    tbl['flags'] = np.asarray(phot['flags'])
    tbl.sort('id')

    tbl.meta = {
        'image': path.name, 'setup': hdr.get('PNSETUP'), 'comb_id': hdr.get('PNCOMBID'),
        'object': values.get('object'), 'filter': values.get('filter'), 'exptime': exptime,
        'airmass': values.get('airmass'), 'mjd': mjd,
        'psf': PSFModel.from_fit(fit).to_meta(scale, nstack),
        'aper_radius_pix': float(settings['aper_radius'] * fwhm), 'removed': removed,
        'settings': settings, 'settings_hash': json_hash(settings),
        'image_hash': image_hash(hdr), 'wcs_method': hdr.get('PNASTMTH', 'none'),
        'wcs_hash': solution_hash(hdr), 'pipnick_version': version
    }
    if 'PNASTMTH' in hdr:
        coordinates = 'celestial coordinates from the astrometric solution.'
    else:
        coordinates = 'no astrometric solution, so the celestial coordinates are masked.'
    log.info(
        f'{path.name}: measured {len(tbl)} sources; PSF FWHM {fwhm:.2f} pixels '
        f'({fwhm * scale:.2f} arcsec) from {nstack} stars; {coordinates}'
    )
    return tbl, fit, sub, mask


def vet_sources(phot, shape, fwhm):
    """
    Remove sources whose PSF fits are not credible.

    Sources are removed if (in this order):

        - their fitted position or flux is not finite (the fit failed);
        - their fitted position is outside the image;
        - their fitted position moved more than one FWHM from where they were
          detected, i.e., onto another source (which is at least one FWHM
          away), typically for spurious detections in the residuals of bright
          sources;
        - their fitted flux is not positive.

    Parameters
    ----------
    phot : :class:`astropy.table.Table`
        The results of the PSF photometry; see
        :func:`~pipnick.photometry.psf_photometry.psf_analysis`.
    shape : tuple
        The shape of the image (rows, columns).
    fwhm : float
        The FWHM of the PSF in pixels.

    Returns
    -------
    phot : :class:`astropy.table.Table`
        The remaining sources.
    removed : dict
        The number of sources removed for each reason: ``failed``,
        ``outside``, ``moved``, and ``nonpositive``.
    """
    x = np.asarray(phot['x_fit'], dtype=float)
    y = np.asarray(phot['y_fit'], dtype=float)
    flux = np.asarray(phot['flux_fit'], dtype=float)
    failed = np.logical_not(np.isfinite(x) & np.isfinite(y) & np.isfinite(flux))
    # NOTE: The comparisons below are only made for finite values
    _x = np.where(failed, 0., x)
    _y = np.where(failed, 0., y)
    outside = np.logical_not(failed) & (
        (_x < -0.5) | (_x > shape[1] - 0.5) | (_y < -0.5) | (_y > shape[0] - 0.5)
    )
    moved = np.logical_not(failed | outside) & (
        np.hypot(_x - np.asarray(phot['x_init']), _y - np.asarray(phot['y_init'])) > fwhm
    )
    nonpositive = np.logical_not(failed | outside | moved) & (np.where(failed, 1., flux) <= 0)
    removed = {
        'failed': int(np.sum(failed)), 'outside': int(np.sum(outside)),
        'moved': int(np.sum(moved)), 'nonpositive': int(np.sum(nonpositive))
    }
    if sum(removed.values()) > 0:
        log.info(
            f'Removed {sum(removed.values())} of {len(phot)} sources: '
            + ', '.join(f'{n} {k}' for k, n in removed.items() if n > 0) + '.'
        )
    return phot[np.logical_not(failed | outside | moved | nonpositive)], removed


def photometry_all(
    table, rdxdir=None, images=None, overwrite=False, thresh=8., fwhm=None, fittype='circ',
    mode='all', stack='mean', group=False, min_sep=1., group_sep=2., aper_radius=2.,
    psf_thresh=15., qa=False, report=None
):
    r"""
    Measure the sources in the reduced science images.

    For each image, the PSF is modeled and the sources are measured by PSF and
    aperture photometry; see
    :func:`~pipnick.photometry.psf_photometry.psf_analysis`.  The object,
    filter, exposure time, and airmass are taken from the metadata table.
    The results are written to ``photometry/<image>_phot.ecsv`` in the
    reduction directory, with these columns:

        - ``id``, ``group_id``, ``group_size``, ``iter_detected``: the source
          and the group of sources fit together; see
          :class:`photutils.psf.IterativePSFPhotometry`;
        - ``x``, ``y``, ``x_err``, ``y_err``: the 0-based pixel coordinates
          (always provided);
        - ``ra``, ``dec``: the ICRS coordinates in degrees, if the image has
          an astrometric solution found by pipnick, and masked otherwise;
        - ``flux_psf``, ``flux_psf_err``: the PSF flux in electrons;
        - ``flux_aper``, ``flux_aper_err``: the background-subtracted flux in
          a circular aperture of radius ``aper_radius`` times the FWHM, masked
          if the aperture includes masked pixels;
        - ``mag_inst``, ``mag_inst_err``: the instrumental magnitude,
          :math:`-2.5 \log_{10}(f / t)` for the PSF flux :math:`f` and the
          exposure time :math:`t`;
        - ``local_bkg``, ``qfit``, ``cfit``, ``reduced_chi2``, ``flags``: the
          local background and the fit-quality measures and flags of
          ``photutils``.

    Sources whose PSF fits are not credible are removed; see
    :func:`vet_sources`.  Values that cannot be measured are masked.  The
    table metadata record the image, its metadata, the PSF model, the
    settings, and the provenance of the image and its astrometric solution.

    When the photometry is rerun, an image is measured again only if it, or
    the settings, have changed (or ``overwrite`` is set).  If only its
    astrometric solution has changed (e.g., the astrometry was run after the
    photometry), only the celestial coordinates are updated.

    A summary of each image is logged and written to
    ``photometry/summary.ecsv``.

    Parameters
    ----------
    table : str, :class:`pathlib.Path`
        The per-setup metadata file used to reduce the images; see
        :func:`~pipnick.pipelines.reduction.reduce_all`.
    rdxdir : str, :class:`pathlib.Path`, optional
        The reduction directory.  If None, the directory saved in the metadata
        file is used.
    images : list, optional
        Only measure these images, given by file name, file name without the
        extension, or combination group (``comb_id``).  If None, all images
        are measured.
    overwrite : bool, optional
        Measure the images again, even if they are up to date.
    thresh : float, optional
        The source detection threshold, as a multiple of the background noise.
    fwhm : float, optional
        The initial guess for the FWHM of the PSF in arcsec.  If None, the
        typical value for the camera is used.
    fittype : str, optional
        The PSF model: ``'circ'`` (circular) or ``'ellip'`` (elliptical)
        Moffat.
    mode : str, optional
        The iterative fitting mode, ``'all'`` or ``'new'``; see
        :class:`photutils.psf.IterativePSFPhotometry`.
    stack : str, optional
        How the stars used to model the PSF are stacked: ``'mean'`` or
        ``'median'``.
    group : bool, optional
        Consolidate the sources fit together into single sources; see
        :func:`~pipnick.photometry.psf_photometry.consolidate_groups`.
    min_sep : float, optional
        The minimum separation of detected sources, in units of the FWHM.
    group_sep : float, optional
        Sources closer than this, in units of the FWHM, are fit together.
    aper_radius : float, optional
        The radius of the photometric apertures, in units of the FWHM.
    psf_thresh : float, optional
        The detection threshold, as a multiple of the background noise, for
        the stars used to model the PSF.
    qa : bool, optional
        Write QA figures for each measured image to ``photometry/qa``: the PSF
        model (``<image>_psf.png``; see
        :func:`~pipnick.qa.psf.plot_psf`) and the measured sources
        (``<image>_sources.png``).
    report : :class:`~pipnick.pkg.report.StepReport`, optional
        Collects the failures and problems of the photometry for its closing
        verdict.  If None, a new one is used.

    Returns
    -------
    list
        The paths to the source tables.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickPathError`
        Raised if the reduction directory does not exist.
    """
    _report = StepReport('PHOTOMETRY') if report is None else report
    dataset = DataSet.from_file(table, rdxdir=rdxdir)
    if not dataset.rdxdir.is_dir():
        raise PipnickPathError(f'{dataset.rdxdir} does not exist; reduce the data first.')
    _report.check_version(dataset.metafile.name, dataset.table.meta.get('pipnick_version'))
    phot_dir = dataset.outdir('photometry')
    reduced = dataset.science_images(selection=images, include_missing=True)
    log.info(f'Photometry of {len(reduced)} images in {dataset.rdxdir}')
    settings = {
        'thresh': float(thresh),
        'fwhm': float(dataset.camera.fwhm_guess if fwhm is None else fwhm),
        'fittype': fittype, 'mode': mode, 'stack': stack, 'group': bool(group),
        'min_sep': float(min_sep), 'group_sep': float(group_sep),
        'aper_radius': float(aper_radius), 'psf_thresh': float(psf_thresh)
    }

    outputs = []
    summary = _summary()
    for k, image in enumerate(reduced):
        progress = f'[{k + 1}/{len(reduced)}]'
        path = image.path
        ofile = image.phot_file
        row = summary.add(image=path.name, status='failed', wcs='none', reason='-')
        if not image.exists:
            _report.attention(
                f'{progress} {path.name} does not exist; reduce the data first.  Skipping.'
            )
            row['status'] = 'missing'
            continue
        hdr = image.header
        row['wcs'] = hdr.get('PNASTMTH', 'none')

        # Check if the existing table is up to date
        if ofile.is_file() and not overwrite:
            existing = SourceCatalog.read(ofile)
            if (
                existing.meta.get('image_hash') == image_hash(hdr)
                and existing.meta.get('settings_hash') == json_hash(settings)
            ):
                row.update(_summary_values(existing))
                _report.check_version(ofile.name, existing.meta.get('pipnick_version'))
                if existing.meta.get('wcs_hash') == solution_hash(hdr):
                    log.info(f'{progress} {ofile.name} is up to date.')
                    row['status'] = 'up to date'
                else:
                    existing.update_coordinates(hdr, solution_hash(hdr))
                    existing.write()
                    log.info(
                        f'{progress} The astrometric solution of {path.name} changed; updated '
                        f'the celestial coordinates in {ofile.name}.'
                    )
                    row['status'] = 'coordinates updated'
                outputs.append(ofile)
                continue

        log.info(f'{progress} Measuring {path.name}.')
        try:
            tbl, fit, sub, mask = measure_image(image, settings)
        except Exception as e:
            log.debug(f'Traceback of the failure of {path.name}:', exc_info=True)
            _report.fail(f'Photometry of {path.name} failed: {type(e).__name__}: {e}')
            row['reason'] = f'{type(e).__name__}: {e}'
            continue
        phot_dir.mkdir(parents=True, exist_ok=True)
        tbl.write(ofile, format='ascii.ecsv', overwrite=True)
        log.info(f'Wrote photometry/{ofile.name}')
        row['status'] = 'measured'
        row.update(_summary_values(tbl))
        outputs.append(ofile)
        if qa:
            _write_qa(path, tbl, fit, sub, mask, phot_dir / 'qa')

    expected = dataset.camera.expected_range('fwhm')
    for row in summary.rows:
        if row['fwhm_arcsec'] is not None:
            _report.check_range(
                f'{row["image"]}: the PSF FWHM', row['fwhm_arcsec'], expected, ' arcsec',
                advice='Check the focus and tracking.'
            )
    if summary.write(phot_dir) is None:
        log.warning('No reduced science images to measure.')
    measured = [row for row in summary.rows if row['fwhm_arcsec'] is not None]
    counts = ', '.join(
        f'{sum(row["status"] == s for row in summary.rows)} {s}'
        for s in ['measured', 'up to date', 'coordinates updated', 'missing', 'failed']
        if any(row['status'] == s for row in summary.rows)
    )
    text = f'{len(summary)} images ({counts})'
    if len(measured) > 0:
        text += f'; median FWHM {np.median([r["fwhm_arcsec"] for r in measured]):.2f} arcsec'
    following = f'pipnick_calibrate -t {dataset.metafile} (for standards)'
    if 'mosaic' in dataset.table.colnames and np.any(dataset.table['mosaic'] != ''):
        following += f', or pipnick_mosaic -t {dataset.metafile}'
    _report.finish(f'{text}; summary in photometry/summary.ecsv.', next_step=following)
    return outputs


def _summary_values(tbl):
    """
    Return the summary values of a source table; see :func:`_summary`.
    """
    return {
        'nsources': len(tbl), 'fwhm_pix': tbl.meta['psf']['fwhm_pix'],
        'fwhm_arcsec': tbl.meta['psf']['fwhm_arcsec'], 'nstack': tbl.meta['psf']['nstack']
    }


def _summary():
    """
    Return the summary of the photometry, written to
    ``photometry/summary.ecsv``.

    The summary gives the status of each image (``measured``, ``up to date``,
    ``coordinates updated``, ``missing``, or ``failed``), the number of
    sources, the FWHM of the PSF (pixels and arcsec), the number of stars used
    to model it, the method of the astrometric solution (``none`` if the
    celestial coordinates are masked), and the reason for any failure.
    """
    return Summary(
        'Photometry',
        ['image', 'status', 'nsources', 'fwhm_pix', 'fwhm_arcsec', 'nstack', 'wcs', 'reason'],
        dtypes={'nsources': int, 'fwhm_pix': float, 'fwhm_arcsec': float, 'nstack': int},
        formats={'fwhm_pix': '.2f', 'fwhm_arcsec': '.2f'},
        logged=['image', 'status', 'nsources', 'fwhm_pix', 'fwhm_arcsec', 'nstack', 'wcs']
    )


def _write_qa(path, tbl, fit, data, mask, qa_dir):
    """
    Write the QA figures for an image: the PSF model and the measured sources.
    """
    qa_dir.mkdir(parents=True, exist_ok=True)
    psf = tbl.meta['psf']
    plot_psf(
        qa_dir / f'{path.stem}_psf.png', fit,
        title=f'{path.name}: PSF FWHM {psf["fwhm_pix"]:.2f} pixels ({psf["fwhm_arcsec"]:.2f} '
              f'arcsec) from {psf["nstack"]} stars'
    )
    order = np.argsort(-tbl['flux_psf'].filled(0.))
    plot_sources(
        data, mask, tbl[order], qa_dir / f'{path.stem}_sources.png',
        title=f'{path.name}: {len(tbl)} sources, numbered by decreasing PSF flux'
    )
