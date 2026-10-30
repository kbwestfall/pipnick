"""
Astrometric calibration of the reduced science images.

The astrometric solution (WCS) of each image is found with as few submissions
to the astrometry.net web service as possible:

    1. Existing solutions are kept, and solutions are restored from the cache
       of earlier astrometry.net results and from solutions provided by the
       user (see :func:`~pipnick.astrometry.header.apply_wcs`), e.g., after the
       images are rewritten by ``pipnick_reduce``.

    2. The images are grouped by their pointing (see
       :func:`~pipnick.astrometry.solve.group_by_pointing`).  In each group,
       one image (the reference) is solved, and the others are solved by
       matching their sources to those of the reference (*bootstrapping*).

    3. A group without a solved image is solved, if possible, from a solved
       image in another group (if their fields overlap), or by matching its
       sources to the Gaia catalog, using the solution predicted from a solved
       image in the same setup and the telescope pointing offset.

    4. An image that still cannot be solved is matched to the Gaia catalog
       with a *wide search*, which allows for an incorrect telescope pointing
       in the metadata table (e.g., a stale pointing in the raw header), with
       a stricter acceptance threshold.  Images with the same incorrect
       pointing then use the pointing error measured for the first.

    5. Only an image that cannot be solved in any of these ways is submitted
       to astrometry.net.

Solutions are written to the headers of the reduced images, with ``PNAST*``
keywords that record how each was found.

All solutions found by the pipeline are in the ICRS.  astrometry.net labels
its solutions as FK5 (J2000), and they are relabeled as ICRS; see
:func:`~pipnick.astrometry.services.icrs_header` for why this is correct.
"""

from astropy.table import Table

from pipnick import log
from pipnick.astrometry.cache import anet_input_hash, invalidate_chains, restore_solution
from pipnick.astrometry.services import GaiaCatalog, astrometry_net_hints
from pipnick.astrometry.solve import (
    AstrometryImage, Session, group_by_pointing, group_order, group_reference, handle_result,
    log_pointing_offset, nsources, read_image, solve_group
)
from pipnick.datamodel.dataset import DataSet
from pipnick.image.detect import brightest, find_sources
from pipnick.pkg.exceptions import PipnickPathError
from pipnick.pkg.paths import qa_figure
from pipnick.pkg.report import StepReport
from pipnick.pkg.settings import AstrometrySettings
from pipnick.pkg.step import Attempt, needs_qa, progress
from pipnick.pkg.summary import Summary
from pipnick.qa.sources import plot_sources


def _detect(data, mask, fwhm, thresh, max_sources):
    """
    Detect the brightest point sources in an image, excluding those at the
    border; see :func:`~pipnick.image.detect.find_sources`.

    The masked pixels are set to 0 before the global background is
    subtracted.  The detected sources are part of the key of the cached
    astrometry.net results (see
    :func:`~pipnick.astrometry.cache.anet_input_hash`), so this detection
    must not change.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The image data.
    mask : :class:`numpy.ndarray`
        Boolean array selecting the pixels to ignore.
    fwhm : float
        The expected FWHM of the sources in pixels.
    thresh : float
        The detection threshold, as a multiple of the standard deviation of
        the background.
    max_sources : int
        The maximum number of sources to return.

    Returns
    -------
    :class:`astropy.table.Table`
        The 0-based pixel coordinates (``x`` and ``y``) and the flux (``flux``)
        of the detected sources, ordered by decreasing flux.
    """
    sources = find_sources(
        data, mask=mask, fwhm=fwhm, thresh=thresh, masked='zero', exclude_border=True
    )
    return brightest(sources, max_sources)['x', 'y', 'flux']


def _load_images(dataset, settings):
    """
    Read the reduced science images and detect their sources.

    A warning is issued for images with fewer than 10 detected sources, which
    are likely to be hard to solve.

    Parameters
    ----------
    dataset : :class:`~pipnick.datamodel.dataset.DataSet`
        The data set.
    settings : :class:`~pipnick.pkg.settings.AstrometrySettings`
        The settings, with the FWHM resolved: the detection settings
        (``max_sources``, ``detect_thresh``, ``fwhm``), and whether images
        will be uploaded to astrometry.net instead of source lists
        (``upload_image``).

    Returns
    -------
    list
        The :class:`~pipnick.astrometry.solve.AstrometryImage` objects, ordered
        by combination group.  Images that do not exist or have no coordinates
        are skipped.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if an image is from a different setup.
    """
    images = []
    reduced = dataset.science_images()
    for k, image in enumerate(reduced):
        label = progress(k, len(reduced))
        path = image.path
        data, mask, hdr = read_image(path)
        pointing = image.pointing
        if pointing is None:
            log.warning(f'{path.name} has no coordinates.  Skipping.')
            continue
        scale = image.nominal_scale
        fwhm_pix = settings.fwhm / scale
        sources = _detect(data, mask, fwhm_pix, settings.detect_thresh, settings.max_sources)
        # Allow for a pointing error of up to the image width
        radius = max(data.shape) * scale / 3600.
        hints = astrometry_net_hints(pointing, radius, scale)
        images.append(
            AstrometryImage(
                path, hdr, data.shape, pointing, image.epoch, scale, fwhm_pix, sources, hints,
                anet_input_hash(data, sources, hints, settings.upload_image),
                comb_id=image.product.comb_id,
                image=image
            )
        )
        log.info(f'{label} {path.name}: detected {len(sources)} sources.')
        if len(sources) < 10:
            log.warning(
                f'Only {len(sources)} sources were detected in {path.name}, which may be hard to '
                'solve; see the troubleshooting advice in the documentation.'
            )
    return images


def _suggestion(image, images):
    """
    Suggest what the user can do to solve an image that could not be solved,
    based on why it could not be solved locally and why it could not be
    submitted; see :func:`~pipnick.astrometry.solve._fail`.
    """
    suggestions = {
        'few_matches': 'check the pointing in the metadata table, or increase --search',
        'few_matches_wide': (
            'check the pointing in the metadata table, or increase --search or --wide_search'
        ),
        'rejected': (
            'the matches may be by chance if the pointing is wrong, so check the pointing in '
            'the metadata table; otherwise, consider --max_rms'
        ),
        'no_reference': (
            'solve another image in its pointing group or setup first'
            if any(i.status == 'solved' for i in images)
            else 'no image in the setup is solved, so a submission is needed'
        ),
        'limit': 'allow a submission (--max_submissions)',
        'waiting': 'rerun to check the pending submission',
        'submission_failed': (
            'try --upload_image or provide a solution; use --resolve and --images to resubmit'
        )
    }
    found = []
    if len(image.sources) < 10:
        found.append(
            'few detected sources: try --detect_thresh, --upload_image, or provide a solution '
            '(fit_wcs and apply_wcs)'
        )
    for reason in [image.reason, image.submit_reason]:
        if reason is not None:
            found.append(suggestions[reason[0]])
    return 'see the log' if len(found) == 0 else '; '.join(found)


def _summarize(images, astro_dir):
    """
    Log a summary of the solution of each image, and write it to
    ``astrometry/summary.ecsv``.

    The summary gives, for each image, its status (``solved``, ``unsolved``,
    ``pending``, or ``skipped`` if not selected), the method, the reference,
    the number of matched sources, the RMS residual (arcsec), and the pointing
    error (arcmin; see :func:`~pipnick.astrometry.solve._pointing_error`) of
    its solution, and, for images that are not solved, the reason and a
    suggestion.
    """
    summary = Summary(
        'Astrometry',
        [
            'image', 'status', 'method', 'reference', 'nmatch', 'rms', 'pointing_error',
            'reason', 'suggestion'
        ],
        dtypes={'nmatch': int, 'rms': float, 'pointing_error': float},
        logged=['image', 'status', 'method', 'reference', 'nmatch', 'rms', 'pointing_error']
    )
    for image in images:
        status = (
            image.status if image.status in ['solved', 'pending', 'skipped'] else 'unsolved'
        )
        hdr = image.header if image.status == 'solved' else {}
        reason = '' if status == 'solved' else ' '.join(
            r[1] for r in [image.reason, image.submit_reason] if r is not None
        )
        suggestion = _suggestion(image, images) if status == 'unsolved' else ''
        summary.add(
            image=image.name, status=status, method=hdr.get('PNASTMTH', '-'),
            reference=hdr.get('PNASTREF', '-'), nmatch=hdr.get('PNASTNPT'),
            rms=hdr.get('PNASTRMS'), pointing_error=hdr.get('PNASTPOE'),
            reason=reason if len(reason) > 0 else '-',
            suggestion=suggestion if len(suggestion) > 0 else '-'
        )
    for row in summary.write(astro_dir):
        if row['status'] == 'unsolved':
            log.warning(
                f'{row["image"]} is not solved: {row["reason"]}  Suggestion: '
                f'{row["suggestion"]}.'
            )


def _write_qa(images, astro_dir, report, solved=None):
    """
    Write the QA figures of the images to ``astrometry/qa/<image>.png``; see
    :func:`~pipnick.qa.sources.plot_sources`.

    A figure is redrawn for each image solved in this run, which shows the
    sources matched by the solution, and drawn for any image without one,
    which shows only the detected sources if the image is unsolved or was
    solved in an earlier run.  A failure is recorded in the report, and does
    not stop the astrometry.

    Parameters
    ----------
    images : list
        The :class:`~pipnick.astrometry.solve.AstrometryImage` objects.
    astro_dir : :class:`pathlib.Path`
        The ``astrometry`` directory.
    report : :class:`~pipnick.pkg.report.StepReport`
        Records any failure.
    solved : set, optional
        The names of the images solved in this run.  If None, no image was
        solved (e.g., for a dry run).
    """
    _solved = set() if solved is None else solved
    nwritten = 0
    for image in images:
        ofile = qa_figure(astro_dir, image.path.stem)
        if not needs_qa(True, image.name in _solved, [ofile]):
            continue
        title = f'{image.name}: {len(image.sources)} detected sources'
        if image.status == 'solved':
            title += f'; {image.method} solution'
        attempt = Attempt(report, f'The QA figure of {image.name}', problem=True)
        with attempt:
            data, mask, _ = read_image(image.path)
            plot_sources(data, mask, image.sources, ofile, matched=image.matched, title=title)
        if not attempt.failed:
            nwritten += 1
    if nwritten > 0:
        log.info(f'Wrote {nwritten} QA figures to astrometry/qa.')


def _select(images, selection):
    """
    Select the images to solve.

    Parameters
    ----------
    images : list
        The :class:`~pipnick.astrometry.solve.AstrometryImage` objects.  Their
        ``selected`` attribute is set.
    selection : list
        The images to select, by file name, file name without the extension,
        or combination group (``comb_id``).

    Raises
    ------
    ValueError
        Raised if an entry does not match any image.
    """
    unknown = []
    for image in images:
        image.selected = False
    for item in selection:
        found = False
        for image in images:
            if (
                str(item) in [image.name, image.path.stem]
                or (str(item).isdigit() and image.comb_id == int(item))
            ):
                image.selected = True
                found = True
        if not found:
            unknown.append(str(item))
    if len(unknown) > 0:
        raise ValueError(f'No reduced science image matches: {", ".join(unknown)}.')


def _report_dry_run(images, groups, astro_dir):
    """
    Report what :func:`astrometry_all` would do, and write the source lists
    that would be submitted to astrometry.net.
    """
    outdir = astro_dir / 'dry_run'
    outdir.mkdir(parents=True, exist_ok=True)
    lines = ['Dry run; nothing is submitted.  Pointing groups:']
    nunsolved = 0
    for i, group in enumerate(groups):
        lines.append(f'    Group {i+1}:')
        for image in group:
            state = (
                f'solved ({image.method})' if image.status == 'solved'
                else 'not selected' if image.status == 'skipped'
                else 'not submitted (earlier submission failed)' if image.submit_failed
                else 'unsolved'
            )
            lines.append(f'        {image.name}: {len(image.sources)} sources; {state}')
            tbl = Table(
                {'x': image.sources['x'], 'y': image.sources['y'], 'flux': image.sources['flux']},
                meta={'hints': image.hints, 'input_hash': image.input_hash}
            )
            tbl.write(outdir / f'{image.path.stem}_sources.ecsv', overwrite=True)
        if group_reference(group) is None:
            nunsolved += 1
            candidates = sorted(group, key=nsources, reverse=True)
            lines.append(f'        Reference if a submission is needed: {candidates[0].name}')
    order = [
        max(g, key=nsources).name for g in sorted(groups, key=group_order)
        if group_reference(g) is None
    ]
    if len(order) > 0:
        lines.append(f'    Order of any submissions: {", ".join(order)}')
    nsolved = len([i for i in images if i.status == 'solved'])
    minimum = 1 if nsolved == 0 and nunsolved > 0 else 0
    lines.append(
        f'{nunsolved} groups without a solved image.  At least {minimum} and, if all '
        f'bootstrapping fails, at most {len([i for i in images if i.status is None])} '
        'astrometry.net submissions would be made.'
    )
    lines.append(f'Source lists written to {outdir}.')
    log.info('\n'.join(lines))


def astrometry_all(
    table, rdxdir=None, solver=None, catalog=None, api_key=None, resolve=False,
    upload_image=False, dry_run=False, max_submissions=None, group_radius=1., min_matches=4,
    max_rms=0.3, search=60., wide_search=10., wide_min_matches=10, images=None, qa=True,
    max_sources=100, detect_thresh=5., fwhm=None, report=None
):
    """
    Find the astrometric solutions of the reduced science images.

    See the module docstring for the strategy, which minimizes the number of
    submissions to astrometry.net.  The solutions are written to the image
    headers; see :func:`~pipnick.astrometry.header.apply_wcs`.  The cached
    astrometry.net results, Gaia queries, and user-provided solutions are kept
    in the ``astrometry`` subdirectory of the reduction directory.

    Parameters
    ----------
    table : str, :class:`pathlib.Path`
        The per-setup metadata file used to reduce the images; see
        :func:`~pipnick.pipelines.reduction.reduce_all`.
    rdxdir : str, :class:`pathlib.Path`, optional
        The reduction directory.  If None, the directory saved in the metadata
        file is used.
    solver : object, optional
        The astrometry.net solver.  If None,
        :class:`~pipnick.astrometry.services.AstrometryNetSolver` is
        constructed when it is first needed.
    catalog : object, optional
        The source catalog.  If None,
        :class:`~pipnick.astrometry.services.GaiaCatalog` is used, with its
        queries cached in the ``astrometry`` directory.
    api_key : str, optional
        The astrometry.net API key, if ``solver`` is None; see
        :class:`~pipnick.astrometry.services.AstrometryNetSolver`.
    resolve : bool, optional
        Ignore the existing and cached solutions, and solve the images again
        (only the selected images, if ``images`` is provided).  This can
        submit jobs to astrometry.net for images that were solved before.
    upload_image : bool, optional
        Upload images to astrometry.net, instead of the lists of their
        detected sources.
    dry_run : bool, optional
        Detect the sources, group the images, and report what would be
        submitted, without solving or submitting anything.
    max_submissions : int, optional
        The maximum number of submissions to astrometry.net.  If None, there
        is no limit.
    group_radius : float, optional
        The maximum distance in arcmin between the pointings of images in a
        group; see :func:`~pipnick.astrometry.solve.group_by_pointing`.
    min_matches : int, optional
        The minimum number of matched sources for a bootstrapped solution.
    max_rms : float, optional
        The maximum RMS residual in arcsec for a bootstrapped solution.
    search : float, optional
        The maximum error in arcsec in the expected offset between matched
        images (or between an image and the catalog), e.g., from the error in
        the telescope pointing.
    wide_search : float, optional
        For images that cannot be solved otherwise, the maximum error in
        arcmin in the expected position of the image when matching to the Gaia
        catalog, to allow for an incorrect pointing in the metadata table.  If
        0, this wide search is not done.
    wide_min_matches : int, optional
        The minimum number of matched sources for a solution found by the wide
        search.
    images : list, optional
        Only solve these images, given by file name, file name without the
        extension, or combination group (``comb_id``).  The other images are
        not changed, but their solutions are used as references.  If None, all
        images are solved.
    qa : bool, optional
        Write the QA figures to ``astrometry/qa``: a figure is redrawn for
        each image solved in this run, and drawn for any image without one;
        see :func:`~pipnick.qa.sources.plot_sources`.
    max_sources : int, optional
        The maximum number of (brightest) sources detected in each image, and
        submitted to astrometry.net.
    detect_thresh : float, optional
        The source detection threshold, as a multiple of the standard
        deviation of the background.
    fwhm : float, optional
        The FWHM of the sources in arcsec, for their detection.  If None, the
        typical value for the camera is used.
    report : :class:`~pipnick.pkg.report.StepReport`, optional
        Collects the problems of the astrometry for its closing verdict: the
        unsolved images, and solutions with a large RMS (see
        :meth:`~pipnick.cameras.cameras.NickelCamera.expected_range`).  If
        None, a new one is used.

    Changing ``max_sources``, ``detect_thresh``, or ``fwhm`` changes the
    sources that would be submitted, so cached astrometry.net results for the
    previous sources are not reused.

    A summary of the solution of each image is logged and written to
    ``astrometry/summary.ecsv``; see :func:`_summarize`.

    Returns
    -------
    list
        The paths to the solved images.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickPathError`
        Raised if the reduction directory does not exist.
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if a setting is invalid; see
        :class:`~pipnick.pkg.settings.AstrometrySettings`.
    """
    _report = StepReport('ASTROMETRY') if report is None else report
    settings = AstrometrySettings(
        max_sources=max_sources, detect_thresh=detect_thresh, fwhm=fwhm,
        upload_image=upload_image, max_submissions=max_submissions, group_radius=group_radius,
        min_matches=min_matches, max_rms=max_rms, search=search, wide_search=wide_search,
        wide_min_matches=wide_min_matches
    )
    dataset = DataSet.from_file(table, rdxdir=rdxdir)
    if not dataset.rdxdir.is_dir():
        raise PipnickPathError(f'{dataset.rdxdir} does not exist; reduce the data first.')
    _report.check_version(dataset.metafile.name, dataset.table.meta.get('pipnick_version'))
    settings = settings.resolve(dataset.camera)
    astro_dir = dataset.outdir('astrometry')
    log.info(f'Astrometric calibration of the images in {dataset.rdxdir}')

    _images = _load_images(dataset, settings)
    if len(_images) == 0:
        _report.attention('There are no reduced science images to calibrate.')
        _report.finish('no images.', next_step=f'pipnick_reduce {dataset.metafile}')
        return []
    if images is not None:
        _select(_images, images)
    by_name = {image.name: image for image in _images}
    for image in _images:
        if not resolve or not image.selected:
            restore_solution(image, by_name, astro_dir)
    invalidate_chains(_images, by_name)
    for image in _images:
        if not image.selected and image.status != 'solved':
            image.status = 'skipped'
        if image.status == 'solved':
            # An existing solution, restored above
            _report.check_version(image.name, image.header.get('PNASTVER'))
    log_pointing_offset(_images)
    groups = group_by_pointing(_images, settings.group_radius)

    if dry_run:
        _report_dry_run(_images, groups, astro_dir)
        if qa:
            _write_qa(_images, astro_dir, _report)
        _report.note('This was a dry run: nothing was submitted to astrometry.net.')
        _report.finish(
            f'dry run of {len(_images)} images.',
            next_step=f'pipnick_astrometry {dataset.metafile}'
        )
        return [image.path for image in _images if image.status == 'solved']

    session = Session(
        solver, GaiaCatalog(astro_dir) if catalog is None else catalog, api_key, astro_dir,
        settings
    )
    session.pending = set(image.name for image in _images if image.status == 'pending')
    previous = set(image.name for image in _images if image.status == 'solved')
    # Solve the groups with a solved image first, then the richest groups
    groups.sort(key=group_order)
    for group in groups:
        solve_group(group, _images, session)

    # Check the pending submissions, and then solve the remaining groups
    for image in _images:
        if image.status != 'pending':
            continue
        status, header, submission_id = session.solver().monitor(image.submission_id)
        handle_result(image, status, header, submission_id, session)
        if image.status == 'pending':
            log.warning(
                f'The astrometry.net submission {image.submission_id} for {image.name} has not '
                'finished; rerun to check it again.'
            )
    groups.sort(key=group_order)
    for group in groups:
        solve_group(group, _images, session)

    if qa:
        solved = {i.name for i in _images if i.status == 'solved' and i.name not in previous}
        _write_qa(_images, astro_dir, _report, solved=solved)
    log_pointing_offset(_images)
    _summarize(_images, astro_dir)
    solved = [image for image in _images if image.status == 'solved']
    methods = {}
    for image in solved:
        methods[image.method] = methods.get(image.method, 0) + 1
    expected = dataset.camera.expected_range('astrometric_rms')
    for image in solved:
        if image.header.get('PNASTRMS') is not None:
            _report.check_range(
                f'{image.name}: the RMS of the astrometric solution',
                float(image.header['PNASTRMS']), expected, ' arcsec',
                advice='Check its QA figure in astrometry/qa, and consider solving it again '
                       '(--resolve).'
            )
    unsolved = [image.name for image in _images if image.status not in ['solved', 'skipped']]
    if len(unsolved) > 0:
        _report.attention(
            f'{len(unsolved)} image(s) are not solved or pending ({", ".join(unsolved)}); see '
            'the reasons and suggestions above.', warn=False
        )
    _report.finish(
        f'{len(solved)} of {len(_images)} images solved ('
        + ', '.join(f'{n} {m}' for m, n in sorted(methods.items()))
        + f'), with {session.nsubmitted} astrometry.net submissions; summary in '
        'astrometry/summary.ecsv.',
        next_step=f'pipnick_photometry {dataset.metafile}'
    )
    return [image.path for image in solved]
