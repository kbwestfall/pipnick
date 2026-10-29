"""
The bookkeeping shared by the loops of the processing steps.

Each step processes its items (images, nights, or mosaics) in its own loop,
and uses these helpers for what every loop does: label the progress
(:func:`progress`), reuse a product that is up to date (:func:`reuse`),
record a failure and continue with the other items (:class:`Attempt`),
decide whether to draw the QA figures (:func:`needs_qa`), and count the
items by status for the closing verdict (:func:`count_statuses`).
"""

from pathlib import Path

from pipnick import log


def progress(k, n):
    """
    Return the progress label of an item, e.g., ``'[3/17]'``.

    Parameters
    ----------
    k : int
        The 0-based index of the item.
    n : int
        The number of items.

    Returns
    -------
    str
        The label.
    """
    return f'[{k + 1}/{n}]'


def reuse(report, name, version, label=None):
    """
    Log that a product is up to date, and check the pipnick version that made
    it; see :meth:`~pipnick.pkg.report.StepReport.check_version`.

    Parameters
    ----------
    report : :class:`~pipnick.pkg.report.StepReport`
        The report of the step.
    name : str
        The name of the product.
    version : str
        The pipnick version recorded in the product, or None.
    label : str, optional
        The progress label; see :func:`progress`.
    """
    log.info(f'{name} is up to date.' if label is None else f'{label} {name} is up to date.')
    report.check_version(name, version)


class Attempt:
    """
    A context manager that processes one item of a step, recording an
    exception as a failure instead of stopping the step.

    The traceback is logged at the debug level, and the failure is recorded
    in the report and in the summary row of the item.  Exceptions that are
    not :class:`Exception` (e.g., :class:`KeyboardInterrupt`) are not caught.

    Parameters
    ----------
    report : :class:`~pipnick.pkg.report.StepReport`
        The report of the step.
    what : str
        What is attempted, e.g., ``'Photometry of M11_V_C096.fits'``; the
        failure is reported as ``'<what> failed: <error>'``.
    row : dict, optional
        The summary row of the item, whose ``reason`` is set to the error.
    then : str, optional
        What the step continues with, e.g., ``'the other images'``, added to
        the message.
    problem : bool, optional
        Record the failure as a problem that needs attention (e.g., for a QA
        figure), rather than a failure of the step.

    Attributes
    ----------
    failed : bool
        True if the item failed.

    Examples
    --------
    .. code-block:: python

        attempt = Attempt(report, f'Photometry of {name}', row=row)
        with attempt:
            tbl = measure_image(image, settings)
        if attempt.failed:
            continue
    """
    def __init__(self, report, what, row=None, then=None, problem=False):
        self.report = report
        self.what = what
        self.row = row
        self.then = then
        self.problem = problem
        self.failed = False

    def __enter__(self):
        return self

    def __exit__(self, kind, error, traceback):
        if error is None or not isinstance(error, Exception):
            return False
        log.debug(f'Traceback of the failure: {self.what}', exc_info=(kind, error, traceback))
        reason = f'{type(error).__name__}: {error}'
        message = f'{self.what} failed: {reason}'
        if self.then is not None:
            message += f'.  Continuing with {self.then}.'
        if self.problem:
            self.report.attention(message)
        else:
            self.report.fail(message)
        if self.row is not None:
            self.row['reason'] = reason
        self.failed = True
        return True


def needs_qa(qa, rebuilt, figures):
    """
    Determine whether the QA figures of a product must be drawn: they are
    redrawn when the product is rebuilt, and drawn when one is missing.

    Parameters
    ----------
    qa : bool
        Whether the step writes QA figures.
    rebuilt : bool
        Whether the product was (re)built in this run.
    figures : list
        The QA figures of the product (:class:`pathlib.Path`).

    Returns
    -------
    bool
        True if the figures must be drawn.
    """
    return qa and (rebuilt or not all(Path(f).is_file() for f in figures))


def count_statuses(rows, statuses, keep_zero=False):
    """
    Count the items of a step by status, for its closing verdict.

    Parameters
    ----------
    rows : list
        The summary rows, each with a ``status``.
    statuses : list
        The statuses, in the order they are listed.
    keep_zero : bool, optional
        List the statuses with no items.

    Returns
    -------
    str
        The counts, e.g., ``'12 measured, 7 up to date'``.
    """
    counts = [(sum(row['status'] == s for row in rows), s) for s in statuses]
    return ', '.join(f'{n} {s}' for n, s in counts if keep_zero or n > 0)
