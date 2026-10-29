"""
A static HTML report of the processing of a data set.

The report collects the status of every step for every science image (see
:func:`~pipnick.report.status.night_status`), the provenance of each image
(see :func:`~pipnick.report.records.image_record`), the summary tables and QA
figures of every step, and the logs, in ``<rdxdir>/report``; see
:func:`write_report`.  It only reads the products, which it never changes or
rebuilds.
"""

import html
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from astropy.table import Table

from pipnick import log
from pipnick.datamodel.dataset import DataSet
from pipnick.pkg.exceptions import PipnickPathError
from pipnick.pkg.paths import step_dir
from pipnick.pkg.report import StepReport
from pipnick.pkg.version import version
from pipnick.report.html import (
    gallery_html, record_html, relative, table_html, write_page
)
from pipnick.report.records import image_record
from pipnick.report.status import night_status


def write_report(table, rdxdir=None, mosaic_dir=None, report=None):
    """
    Write a static HTML report of the processing of a data set.

    The report is written to the ``report`` subdirectory of the reduction
    directory:

        - ``index.html``: the status of every step for every science image,
          colored by state, with the images that need attention, the next
          command to run, links to the pages of the steps, and the logs;
        - ``images/<image>.html``: the provenance of each science image (see
          :func:`~pipnick.report.records.image_record`), with its QA figures;
        - ``reduce.html``, ``astrometry.html``, ``photometry.html``,
          ``calibrate.html``, and ``mosaic.html``: the summary table and the QA
          figures of each step.  The calibration page also lists the bands
          whose extinction was not measured.

    Every link is relative, so the report can be opened from a copy of the
    whole reduction directory.  Products and figures that do not exist are
    noted as such; nothing is rebuilt.

    Parameters
    ----------
    table : str, :class:`pathlib.Path`
        A per-setup metadata file.
    rdxdir : str, :class:`pathlib.Path`, optional
        The reduction directory.  If None, the directory saved in the metadata
        file is used.
    mosaic_dir : str, :class:`pathlib.Path`, optional
        The directory with the mosaics.  If None, the ``mosaic`` subdirectory
        of the reduction directory is used.
    report : :class:`~pipnick.pkg.report.StepReport`, optional
        Collects the problems for the closing verdict.  If None, a new one is
        used.

    Returns
    -------
    :class:`pathlib.Path`
        The main page of the report, ``index.html``.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickPathError`
        Raised if the reduction directory does not exist.
    """
    _report = StepReport('REPORT') if report is None else report
    dataset = DataSet.from_file(table, rdxdir=rdxdir, read_only=True)
    rdx = dataset.rdxdir.absolute()
    if not rdx.is_dir():
        raise PipnickPathError(f'{rdx} does not exist; reduce the data first.')
    _mosaic_dir = step_dir(rdx, 'mosaic') if mosaic_dir is None else Path(mosaic_dir).absolute()
    root = step_dir(rdx, 'report')
    subtitle = (
        f'Setup {dataset.tag} ({dataset.camera.name}); metadata {dataset.metafile.name}; '
        f'reduction directory {rdx}; pipnick {version}; written '
        f'{datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")} UT'
    )
    log.info(f'Writing the report of {rdx} to {root}')

    status, next_step = night_status(table, rdxdir=rdx, mosaic_dir=_mosaic_dir)
    images = dataset.science_images(include_missing=True)
    by_name = {image.name: image for image in images}
    image_dir = root / 'images'
    for image in images:
        record = image_record(image, by_name=by_name, mosaic_dir=_mosaic_dir)
        page = image_dir / f'{image.path.stem}.html'
        write_page(page, image.name, record_html(record, page.parent), root, subtitle=subtitle)
    links = {
        'image': {image.name: f'images/{image.path.stem}.html' for image in images}
    }
    _write_index(root, status, next_step, links, rdx, subtitle)

    steps = [
        ('reduce', 'Reduction', step_dir(rdx, 'reduce')),
        ('astrometry', 'Astrometry', step_dir(rdx, 'astrometry')),
        ('photometry', 'Photometry', step_dir(rdx, 'photometry')),
        ('calibrate', 'Calibration', step_dir(rdx, 'calibrate')),
        ('mosaic', 'Mosaics', _mosaic_dir)
    ]
    missing = []
    for name, title, outdir in steps:
        if not (outdir / 'summary.ecsv').is_file():
            missing.append(title.lower())
        _write_step(root / f'{name}.html', title, outdir, links, root, subtitle)
    if len(missing) > 0:
        _report.note(f'Steps without a summary (not yet run): {", ".join(missing)}.')
    attention = int(np.sum(np.asarray(status['attention'].filled('')) != ''))
    if attention > 0:
        _report.attention(
            f'{attention} image(s) need attention; see the attention column of the report.',
            warn=False
        )
    _report.finish(
        f'{len(images)} images; open {root / "index.html"}', next_step=next_step
    )
    return root / 'index.html'


def _write_index(root, status, next_step, links, rdxdir, subtitle):
    """
    Write the main page of the report; see :func:`write_report`.
    """
    states = ['reduced', 'wcs', 'phot', 'calibrated', 'mosaic']
    parts = [
        f'<p><b>Next:</b> {html.escape(next_step)}</p>',
        '<h2>Status of the science images</h2>',
        '<p>Each step is <span class="good">ok</span>, <span class="warn">stale</span> (run '
        'it again), <span class="muted">missing</span>, <span class="bad">failed</span>, or '
        'n/a.  Follow an image to its provenance and QA figures.</p>',
        table_html(status, links=links, status_columns=states)
    ]
    attention = status[np.asarray(status['attention'].filled('')) != '']
    if len(attention) > 0:
        items = ''.join(
            f'<li><a href="{html.escape(links["image"][str(row["image"])])}">'
            f'{html.escape(str(row["image"]))}</a>: {html.escape(str(row["attention"]))}</li>'
            for row in attention
        )
        parts.append(f'<h2>Needs attention</h2>\n<ul>{items}</ul>')
    logs = sorted(step_dir(rdxdir, 'log').glob('*.log'), reverse=True)
    if len(logs) > 0:
        items = ''.join(
            f'<li><a href="{html.escape(relative(f, root))}">{html.escape(f.name)}</a></li>'
            for f in logs
        )
        parts.append(f'<h2>Logs (newest first)</h2>\n<ul>{items}</ul>')
    write_page(root / 'index.html', 'pipnick report', '\n'.join(parts), root, subtitle=subtitle)


def _write_step(path, title, outdir, links, root, subtitle):
    """
    Write the page of a step: its summary table and QA figures; see
    :func:`write_report`.
    """
    summary_file = outdir / 'summary.ecsv'
    if not summary_file.is_file():
        body = f'<p class="muted">There is no summary in {html.escape(str(outdir))}.</p>'
        write_page(path, title, body, root, subtitle=subtitle)
        return
    summary = Table.read(summary_file)
    parts = [
        f'<p>Summary table: <a href="{html.escape(relative(summary_file, root))}">'
        f'{html.escape(summary_file.name)}</a></p>'
    ]
    if 'extinction' in summary.colnames:
        fallback = summary[np.asarray(summary['extinction']) == 'fallback']
        if len(fallback) > 0:
            items = ''.join(
                f'<li>{html.escape(str(r["night"]))} {html.escape(str(r["band"]))}: '
                f'{html.escape(str(r["reason"]))}</li>'
                for r in fallback
            )
            parts.append(
                '<div class="banner"><b>Extinction not measured:</b> the nominal extinction '
                f'was used for these bands, without its systematic error.<ul>{items}</ul></div>'
            )
    columns = [c for c in ['status', 'action'] if c in summary.colnames]
    parts.append(table_html(summary, links=links, status_columns=columns))
    figures = sorted((outdir / 'qa').glob('*.png'))
    if len(figures) > 0:
        parts.append('<h2>QA figures</h2>')
        parts.append(gallery_html(figures, root))
    write_page(path, title, '\n'.join(parts), root, subtitle=subtitle)
