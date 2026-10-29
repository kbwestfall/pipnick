"""
Rendering of the pages of the report of a data set; see
:func:`~pipnick.pipelines.report.write_report`.

The pages are static HTML, with an embedded style sheet and no scripts.
Every link (to the other pages, the QA figures, the summary tables, and the
logs) is relative, so the report works as long as the reduction directory is
kept together.
"""

import html
from pathlib import Path


def relative(path, page_dir):
    """
    Return the link from a page to a file.

    Parameters
    ----------
    path : str, :class:`pathlib.Path`
        The file.
    page_dir : str, :class:`pathlib.Path`
        The directory of the page.

    Returns
    -------
    str
        The relative path, with ``/`` separators.
    """
    return Path(path).absolute().relative_to(Path(page_dir).absolute(), walk_up=True).as_posix()


def status_class(value):
    """
    Return the style class of a status, used to color it: ``good``
    (e.g., ``ok``, ``up to date``, ``solved``), ``warn`` (e.g., ``stale``,
    ``pending``), ``bad`` (``failed``, ``unsolved``), ``muted`` (``missing``,
    ``n/a``), or an empty string for other values.

    Parameters
    ----------
    value : str
        The status, from the status table (see
        :func:`~pipnick.report.status.night_status`) or the summary of a step.

    Returns
    -------
    str
        The style class.
    """
    classes = {
        'good': [
            'ok', 'up to date', 'solved', 'built', 'measured', 'coordinates updated',
            'constructed', 'fit'
        ],
        'warn': ['stale', 'pending', 'skipped', 'edited (kept)'],
        'bad': ['failed', 'unsolved'],
        'muted': ['missing', 'n/a']
    }
    for name, values in classes.items():
        if str(value) in values:
            return name
    return ''


def table_html(tbl, links=None, status_columns=()):
    """
    Render a table.

    Parameters
    ----------
    tbl : :class:`astropy.table.Table`
        The table; masked values are shown as ``--``, and the column formats
        are used.
    links : dict, optional
        Links for the cells of a column: the column name, and a dictionary
        with the link of each cell value.
    status_columns : tuple, optional
        The columns whose values are colored by :func:`status_class`.

    Returns
    -------
    str
        The table.
    """
    _links = {} if links is None else links
    head = ''.join(f'<th>{html.escape(name)}</th>' for name in tbl.colnames)
    # The cells formatted by astropy, with the column formats and masked
    # values shown as --
    text = {
        name: [
            t.strip()
            for t in tbl[name].pformat(show_name=False, show_unit=False, max_lines=-1)
        ]
        for name in tbl.colnames
    }
    rows = []
    for i in range(len(tbl)):
        cells = []
        for name in tbl.colnames:
            content = html.escape(text[name][i])
            link = _links.get(name, {}).get(text[name][i])
            if link is not None:
                content = f'<a href="{html.escape(link)}">{content}</a>'
            css = status_class(text[name][i]) if name in status_columns else ''
            attr = f' class="{css}"' if len(css) > 0 else ''
            cells.append(f'<td{attr}>{content}</td>')
        rows.append(f'<tr>{"".join(cells)}</tr>')
    return f'<table>\n<tr>{head}</tr>\n' + '\n'.join(rows) + '\n</table>'


def figure_html(caption, path, page_dir):
    """
    Render a QA figure, or a note that it does not exist.

    Parameters
    ----------
    caption : str
        The caption.
    path : :class:`pathlib.Path`
        The figure.
    page_dir : :class:`pathlib.Path`
        The directory of the page.

    Returns
    -------
    str
        The figure, linked to the full-size image.
    """
    if not path.is_file():
        return (
            f'<figure><figcaption>{html.escape(caption)}</figcaption>'
            f'<p class="muted">{html.escape(path.name)} has not been written.</p></figure>'
        )
    link = html.escape(relative(path, page_dir))
    return (
        f'<figure><figcaption>{html.escape(caption)}</figcaption>'
        f'<a href="{link}"><img src="{link}" alt="{html.escape(path.name)}" loading="lazy">'
        '</a></figure>'
    )


def gallery_html(paths, page_dir):
    """
    Render a set of QA figures, each captioned with its file name.

    Parameters
    ----------
    paths : list
        The figures (:class:`pathlib.Path`).
    page_dir : :class:`pathlib.Path`
        The directory of the page.

    Returns
    -------
    str
        The figures.
    """
    figures = ''.join(figure_html(p.stem, p, page_dir) for p in paths)
    return f'<div class="gallery">{figures}</div>'


def record_html(record, page_dir):
    """
    Render a record of an image or a source; see
    :class:`~pipnick.report.records.Record`.

    Parameters
    ----------
    record : :class:`~pipnick.report.records.Record`
        The record.
    page_dir : :class:`pathlib.Path`
        The directory of the page.

    Returns
    -------
    str
        The sections, each as a table of labels and values, and the QA
        figures.
    """
    parts = []
    for title, entries in record.sections.items():
        rows = []
        for label, value in entries:
            css = status_class(value) if label == 'status' else ''
            attr = f' class="{css}"' if len(css) > 0 else ''
            rows.append(
                f'<tr><th>{html.escape(label)}</th><td{attr}>{html.escape(value)}</td></tr>'
            )
        parts.append(
            f'<h2>{html.escape(title)}</h2>\n<table class="record">\n' + '\n'.join(rows)
            + '\n</table>'
        )
    if len(record.figures) > 0:
        figures = ''.join(figure_html(c, p, page_dir) for c, p in record.figures)
        parts.append(f'<h2>QA figures</h2>\n<div class="figures">{figures}</div>')
    return '\n'.join(parts)


def write_page(path, title, body, root, subtitle=None):
    """
    Write a page of the report.

    Parameters
    ----------
    path : :class:`pathlib.Path`
        The page.  Its directory is created if necessary.
    title : str
        The title.
    body : str
        The HTML content of the page.
    root : :class:`pathlib.Path`
        The directory of the report, with ``index.html``, used for the
        navigation links.
    subtitle : str, optional
        A line below the title, e.g., the data set.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    nav = ' | '.join(
        f'<a href="{relative(root / f"{name}.html", path.parent)}">{label}</a>'
        for name, label in [
            ('index', 'Overview'), ('reduce', 'Reduction'), ('astrometry', 'Astrometry'),
            ('photometry', 'Photometry'), ('calibrate', 'Calibration'), ('mosaic', 'Mosaics')
        ]
    )
    sub = '' if subtitle is None else f'<p class="subtitle">{html.escape(subtitle)}</p>'
    style = (
        'body { font-family: sans-serif; margin: 1em 2em; color: #222; }\n'
        'nav { margin-bottom: 1em; }\n'
        'h1 { margin-bottom: 0.2em; }\n'
        '.subtitle { color: #555; margin-top: 0; }\n'
        'table { border-collapse: collapse; margin: 0.5em 0 1.5em 0; font-size: 0.9em; }\n'
        'th, td { border: 1px solid #ccc; padding: 0.2em 0.5em; text-align: left; }\n'
        'th { background: #f0f0f0; }\n'
        'table.record th { width: 14em; }\n'
        '.good { background: #d8f0d8; }\n'
        '.warn { background: #fbe7c2; }\n'
        '.bad { background: #f6caca; }\n'
        '.muted { color: #888; }\n'
        '.banner { background: #f6caca; border: 1px solid #c33; padding: 0.5em 1em; }\n'
        '.figures figure, .gallery figure '
        '{ display: inline-block; margin: 0.5em; vertical-align: top; }\n'
        '.figures img { max-width: 900px; }\n'
        '.gallery img { max-width: 420px; }\n'
        'figcaption { font-weight: bold; margin-bottom: 0.3em; }\n'
    )
    path.write_text(
        '<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        f'<title>{html.escape(title)}</title>\n<style>\n{style}</style>\n</head>\n<body>\n'
        f'<nav>{nav}</nav>\n<h1>{html.escape(title)}</h1>\n{sub}\n{body}\n</body>\n</html>\n'
    )

