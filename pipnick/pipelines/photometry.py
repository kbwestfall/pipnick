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

from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.table import MaskedColumn, Table

from pipnick import log
from pipnick.astrometry.fit import to_icrs, wcs_from_header
from pipnick.astrometry.qa import plot_solution
from pipnick.photometry.aperture_photometry import aperture_analysis
from pipnick.photometry.fit import psf_plot
from pipnick.photometry.psf_photometry import consolidate_groups, psf_analysis
from pipnick.pipelines.reduction import load_setup, plan_products, table_values
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickPathError
from pipnick.pkg.paths import photometry_table, step_dir
from pipnick.pkg.provenance import image_hash, json_hash, solution_hash


def celestial_coordinates(header, x, y):
    """
    Compute the celestial coordinates of sources, if the image has an
    astrometric solution found by pipnick.

    Parameters
    ----------
    header : :class:`astropy.io.fits.Header`
        The header of the reduced image.
    x, y : :class:`numpy.ndarray`
        The 0-based pixel coordinates of the sources.

    Returns
    -------
    ra, dec : :class:`astropy.table.MaskedColumn`
        The ICRS right ascension and declination in degrees, masked if the
        image has no solution.
    """
    n = len(x)
    if 'PNASTMTH' not in header or n == 0:
        masked = np.ones(n, dtype=bool)
        return (
            MaskedColumn(np.zeros(n), name='ra', unit='deg', mask=masked),
            MaskedColumn(np.zeros(n), name='dec', unit='deg', mask=masked)
        )
    coords = to_icrs(wcs_from_header(header).pixel_to_world(x, y))
    unmasked = np.zeros(n, dtype=bool)
    return (
        MaskedColumn(coords.ra.deg, name='ra', unit='deg', mask=unmasked),
        MaskedColumn(coords.dec.deg, name='dec', unit='deg', mask=unmasked)
    )


def _masked(values):
    """
    Return a masked array with the non-finite values (e.g., from failed fits)
    masked and replaced by 0.
    """
    _values = np.asarray(values, dtype=float)
    bad = np.logical_not(np.isfinite(_values))
    return np.ma.MaskedArray(np.where(bad, 0., _values), mask=bad)


def measure_image(path, values, mjd, settings, scale):
    """
    Measure the sources in a reduced image.

    Parameters
    ----------
    path : :class:`pathlib.Path`
        The reduced image.
    values : dict
        The metadata table values for the image (object, filter, exposure
        time, airmass); see :func:`~pipnick.pipelines.reduction.table_values`.
    mjd : float
        The MJD of the (first) observation.
    settings : dict
        The photometry settings; see :func:`photometry_all`.
    scale : float
        The pixel scale in arcsec per pixel.

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
    ccd = CCDData.read(path)
    hdr = fits.getheader(path)
    mask = np.zeros(ccd.data.shape, dtype=bool) if ccd.mask is None else ccd.mask
    # The pixel errors propagated by the reduction (see
    # pipnick.pipelines.reduction.process_frame)
    if ccd.uncertainty is None:
        log.warning(f'{path.name} has no uncertainty; the flux errors are not weighted.')
        error = None
    else:
        error = np.asarray(ccd.uncertainty.represent_as(StdDevUncertainty).array, dtype=float)
        # Pixels without a valid error are ignored
        bad = np.logical_not(np.isfinite(error)) | (error <= 0)
        mask = mask | bad
        error = np.where(bad, 1., error)
    data = np.ma.MaskedArray(np.asarray(ccd.data, dtype=float), mask=mask)

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
        'psf': {
            'fittype': settings['fittype'], 'par': [float(p) for p in fit.par],
            'fwhm_pix': float(fwhm), 'fwhm_arcsec': float(fwhm * scale), 'nstack': nstack
        },
        'aper_radius_pix': float(settings['aper_radius'] * fwhm), 'removed': removed,
        'settings': settings, 'settings_hash': json_hash(settings),
        'image_hash': image_hash(hdr), 'wcs_method': hdr.get('PNASTMTH', 'none'),
        'wcs_hash': solution_hash(hdr)
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


def select_products(products, selection):
    """
    Select the products to process, by file name, file name without the
    extension, or combination group (``comb_id``).

    Parameters
    ----------
    products : list
        The products; see :func:`~pipnick.pipelines.reduction.plan_products`.
    selection : list
        The file names, file names without the extension, or combination
        groups of the products to select.

    Returns
    -------
    list
        The selected products, in their original order.

    Raises
    ------
    :exc:`ValueError`
        Raised if an entry does not match any product.
    """
    selected = []
    unknown = []
    for item in selection:
        found = [
            p for p in products
            if str(item) in [p['file'].name, p['file'].stem]
            or (str(item).isdigit() and p['comb_id'] == int(item))
        ]
        if len(found) == 0:
            unknown.append(str(item))
        selected += [p for p in found if p not in selected]
    if len(unknown) > 0:
        raise ValueError(f'No reduced science image matches: {", ".join(unknown)}.')
    return [p for p in products if p in selected]


def photometry_all(
    table, rdxdir=None, images=None, overwrite=False, thresh=8., fwhm=None, fittype='circ',
    mode='all', stack='mean', group=False, min_sep=1., group_sep=2., aper_radius=2.,
    psf_thresh=15., qa=False
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
        :func:`~pipnick.photometry.fit.psf_plot`) and the measured sources
        (``<image>_sources.png``).

    Returns
    -------
    list
        The paths to the source tables.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickPathError`
        Raised if the reduction directory does not exist.
    """
    metafile, metadata, _rdxdir, tag, camera = load_setup(table, rdxdir=rdxdir)
    if not _rdxdir.is_dir():
        raise PipnickPathError(f'{_rdxdir} does not exist; reduce the data first.')
    phot_dir = step_dir(_rdxdir, 'photometry')
    log.info(f'Photometry of the images in {_rdxdir}')

    products = [p for p in plan_products(metadata, _rdxdir) if p['frametype'] == 'science']
    if images is not None:
        products = select_products(products, images)
    settings = {
        'thresh': float(thresh), 'fwhm': float(camera.fwhm_guess if fwhm is None else fwhm),
        'fittype': fittype, 'mode': mode, 'stack': stack, 'group': bool(group),
        'min_sep': float(min_sep), 'group_sep': float(group_sep),
        'aper_radius': float(aper_radius), 'psf_thresh': float(psf_thresh)
    }

    outputs = []
    summary = []
    for product in products:
        path = product['file']
        ofile = photometry_table(path)
        row = {
            'image': path.name, 'status': 'failed', 'nsources': -1, 'fwhm_pix': -1.,
            'fwhm_arcsec': -1., 'nstack': -1, 'wcs': 'none', 'reason': '-'
        }
        summary.append(row)
        if not path.is_file():
            log.warning(f'{path.name} does not exist; reduce the data first.  Skipping.')
            row['status'] = 'missing'
            continue
        hdr = fits.getheader(path)
        if hdr.get('PNSETUP') != tag:
            raise PipnickDataModelError(f'{path.name} is not from setup {tag}.')
        row['wcs'] = hdr.get('PNASTMTH', 'none')
        i = product['indx'][0]
        scale = camera.plate_scale * int(metadata['xbin'][i])

        # Check if the existing table is up to date
        if ofile.is_file() and not overwrite:
            existing = Table.read(ofile, format='ascii.ecsv')
            if (
                existing.meta.get('image_hash') == image_hash(hdr)
                and existing.meta.get('settings_hash') == json_hash(settings)
            ):
                row.update(_summary_values(existing))
                if existing.meta.get('wcs_hash') == solution_hash(hdr):
                    log.info(f'{ofile.name} is up to date.')
                    row['status'] = 'up to date'
                else:
                    existing['ra'], existing['dec'] = celestial_coordinates(
                        hdr, np.asarray(existing['x']), np.asarray(existing['y'])
                    )
                    existing.meta['wcs_method'] = hdr.get('PNASTMTH', 'none')
                    existing.meta['wcs_hash'] = solution_hash(hdr)
                    existing.write(ofile, format='ascii.ecsv', overwrite=True)
                    log.info(
                        f'The astrometric solution of {path.name} changed; updated the '
                        f'celestial coordinates in {ofile.name}.'
                    )
                    row['status'] = 'coordinates updated'
                outputs.append(ofile)
                continue

        mjd = None if metadata['mjd'][i] is np.ma.masked else float(metadata['mjd'][i])
        try:
            tbl, fit, sub, mask = measure_image(
                path, table_values(metadata, product['indx'], 'science'), mjd, settings, scale
            )
        except ValueError as e:
            log.warning(f'Photometry of {path.name} failed: {e}')
            row['reason'] = str(e)
            continue
        phot_dir.mkdir(parents=True, exist_ok=True)
        tbl.write(ofile, format='ascii.ecsv', overwrite=True)
        log.info(f'Wrote {ofile}')
        row['status'] = 'measured'
        row.update(_summary_values(tbl))
        outputs.append(ofile)
        if qa:
            _write_qa(path, tbl, fit, sub, mask, phot_dir / 'qa')

    _summarize(summary, phot_dir)
    return outputs


def _summary_values(tbl):
    """
    Return the summary values of a source table; see :func:`_summarize`.
    """
    return {
        'nsources': len(tbl), 'fwhm_pix': tbl.meta['psf']['fwhm_pix'],
        'fwhm_arcsec': tbl.meta['psf']['fwhm_arcsec'], 'nstack': tbl.meta['psf']['nstack']
    }


def _summarize(summary, phot_dir):
    """
    Log the summary of the photometry of each image, and write it to
    ``photometry/summary.ecsv``.

    The summary gives the status of each image (``measured``, ``up to date``,
    ``coordinates updated``, ``missing``, or ``failed``), the number of
    sources, the FWHM of the PSF (pixels and arcsec), the number of stars used
    to model it, the method of the astrometric solution (``none`` if the
    celestial coordinates are masked), and the reason for any failure.
    """
    if len(summary) == 0:
        log.warning('No reduced science images to measure.')
        return
    tbl = Table(rows=summary)
    # Missing values are masked, not set to an invalid number
    for name in ['nsources', 'fwhm_pix', 'fwhm_arcsec', 'nstack']:
        tbl[name] = np.ma.MaskedArray(tbl[name], mask=np.asarray(tbl[name]) < 0)
    for name in ['fwhm_pix', 'fwhm_arcsec']:
        tbl[name].info.format = '.2f'
    phot_dir.mkdir(parents=True, exist_ok=True)
    tbl.write(phot_dir / 'summary.ecsv', format='ascii.ecsv', overwrite=True)
    columns = ['image', 'status', 'nsources', 'fwhm_pix', 'fwhm_arcsec', 'nstack', 'wcs']
    log.info(
        'Photometry summary (also in photometry/summary.ecsv):\n'
        + '\n'.join(tbl[columns].pformat(max_lines=-1, max_width=-1))
    )
    for row in tbl:
        if row['status'] == 'failed':
            log.warning(f'Photometry of {row["image"]} failed: {row["reason"]}')


def _write_qa(path, tbl, fit, data, mask, qa_dir):
    """
    Write the QA figures for an image: the PSF model and the measured sources.
    """
    qa_dir.mkdir(parents=True, exist_ok=True)
    psf = tbl.meta['psf']
    psf_plot(
        qa_dir / f'{path.stem}_psf.png', fit,
        title=f'{path.name}: PSF FWHM {psf["fwhm_pix"]:.2f} pixels ({psf["fwhm_arcsec"]:.2f} '
              f'arcsec) from {psf["nstack"]} stars'
    )
    order = np.argsort(-tbl['flux_psf'].filled(0.))
    plot_solution(
        data, mask, tbl[order], qa_dir / f'{path.stem}_sources.png',
        title=f'{path.name}: {len(tbl)} sources, numbered by decreasing PSF flux'
    )
