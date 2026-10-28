"""
The closing verdict of a pipnick step.

Each step collects what it finds as it runs in a :class:`StepReport`:

    - **failures**: products that could not be constructed;
    - **problems**: results that need attention, e.g., a measurement outside
      the expected range of the camera, or a missing calibration;
    - **notes**: information that does not need action;
    - the pipnick **versions** of the existing products that the step
      reuses, so that products made by another version can be reported.

At the end, :meth:`StepReport.finish` logs the verdict: ``<STEP> OK``,
``<STEP> NEEDS ATTENTION`` (with the problems), or ``<STEP> FAILED`` (with the
failures), and the next command to run.  A failed step should end with a
non-zero exit status; see :attr:`StepReport.status`.
"""

from pipnick import log
from pipnick.pkg.version import version


class StepReport:
    """
    The verdict of a pipnick step.

    Parameters
    ----------
    step : str
        The name of the step, used in the verdict (e.g., ``'REDUCTION'``).
    """
    def __init__(self, step):
        self.step = step
        self.failures = []
        self.problems = []
        self.notes = []
        # The existing products and the pipnick version that made them
        self.versions = {}

    def fail(self, message):
        """
        Record (and log) a failure: a product that could not be constructed.

        Parameters
        ----------
        message : str
            The failure.
        """
        log.warning(message)
        self.failures.append(message)

    def attention(self, message, warn=True):
        """
        Record a problem that needs attention.

        Parameters
        ----------
        message : str
            The problem.
        warn : bool, optional
            Also log the problem as a warning.  Set to False if it has
            already been logged.
        """
        if warn:
            log.warning(message)
        self.problems.append(message)

    def note(self, message):
        """
        Record a note for the verdict, which does not need action.

        Parameters
        ----------
        message : str
            The note.
        """
        self.notes.append(message)

    def check_version(self, name, product_version):
        """
        Record the pipnick version of an existing product that is reused.

        Products made by a version other than the installed one are reported
        by :meth:`finish`.

        Parameters
        ----------
        name : str
            The name of the product.
        product_version : str
            The version recorded in the product, or None if it has none
            (i.e., it was made by an earlier version of pipnick).
        """
        self.versions[name] = product_version

    def check_range(self, label, value, limits, unit='', advice=None):
        """
        Check that a measurement is in its expected range, and record a
        problem if it is not.

        Parameters
        ----------
        label : str
            What was measured, e.g., ``'M11_V_C014: PSF FWHM'``.
        value : float
            The measurement.
        limits : tuple
            The lower and upper limits; either can be None.  If ``limits`` is
            None, the measurement is not checked.
        unit : str, optional
            The unit, used in the message.
        advice : str, optional
            What to do if the measurement is out of range.

        Returns
        -------
        bool
            False if the measurement is outside the range; True otherwise.
        """
        if limits is None or value is None:
            return True
        lo, hi = limits
        if (lo is None or value >= lo) and (hi is None or value <= hi):
            return True
        expected = _range_text(lo, hi, unit)
        message = f'{label} is {value:.2f}{unit}, outside the expected range ({expected}).'
        if advice is not None:
            message += f'  {advice}'
        self.attention(message)
        return False

    @property
    def mismatched(self):
        """The existing products made by a version other than the installed one."""
        return {k: v for k, v in self.versions.items() if v != version}

    @property
    def status(self):
        """``'FAILED'``, ``'NEEDS ATTENTION'``, or ``'OK'``."""
        if len(self.failures) > 0:
            return 'FAILED'
        if len(self.problems) > 0:
            return 'NEEDS ATTENTION'
        return 'OK'

    @property
    def exit_status(self):
        """The exit status of a script: 1 if the step failed, 0 otherwise."""
        return 1 if self.status == 'FAILED' else 0

    def version_warning(self):
        """
        Warn if any existing product was made by a version other than the
        installed one, and return the message (or None).
        """
        mismatched = self.mismatched
        if len(mismatched) == 0:
            return None
        names = sorted(mismatched)
        listed = ', '.join(names[:5]) + (f', and {len(names) - 5} more' if len(names) > 5 else '')
        made = sorted(
            {'an earlier version (none recorded)' if v is None else v for v in mismatched.values()}
        )
        message = (
            f'{len(names)} existing product(s) were made with a different version of pipnick '
            f'({", ".join(made)}) than the one installed ({version}): {listed}.  They are still '
            'used, but rerunning all of the processing steps with --overwrite is recommended.'
        )
        log.warning(message)
        return message

    def finish(self, summary, next_step=None):
        """
        Log the verdict of the step.

        Parameters
        ----------
        summary : str
            A summary of what the step did, e.g., ``'17 products (3 up to
            date, 14 built)'``.
        next_step : str, optional
            The next command to run, given unless the step failed.

        Returns
        -------
        str
            The status; see :attr:`status`.
        """
        version_message = self.version_warning()
        status = self.status
        lines = [f'{self.step} {status}: {summary}']
        if status == 'FAILED':
            lines += [f'    FAILED: {m}' for m in self.failures]
        lines += [f'    ATTENTION: {m}' for m in self.problems]
        lines += [f'    NOTE: {m}' for m in self.notes]
        if version_message is not None:
            lines.append(f'    NOTE: {version_message}')
        if status != 'FAILED' and next_step is not None:
            lines.append(f'Next: {next_step}')
        message = '\n'.join(lines)
        if status == 'OK':
            log.info(message)
        else:
            log.warning(message)
        return status


def _range_text(lo, hi, unit):
    """
    Return the description of a range, e.g., ``'0.80-3.50 arcsec'`` or ``'at
    most 0.50 arcsec'``.
    """
    if lo is None:
        return f'at most {hi:.2f}{unit}'
    if hi is None:
        return f'at least {lo:.2f}{unit}'
    return f'{lo:.2f}-{hi:.2f}{unit}'
