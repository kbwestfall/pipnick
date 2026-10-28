"""
Implements base classes for use with scripts.

.. note::

    Adapted from ``pypeit/scripts/scriptbase.py`` in `PypeIt
    <https://github.com/pypeit/PypeIt>`__ (BSD 3-Clause); see
    ``licenses/README.rst``.
"""

import argparse
import datetime
import sys
from functools import partial
from pathlib import Path

from .. import log
from ..pkg.exceptions import PipnickPathError
from ..pkg.version import version


class ScriptBase:
    """
    Provides a base class for all scripts.
    """
    @classmethod
    def entry_point(cls):
        """
        Defines the main script entry point.

        The exit status is the value returned by :meth:`main`, if it is not
        None or 0, e.g., 1 if the step failed.
        """
        status = cls.main(cls.parse_args())
        if status is not None and status != 0:
            sys.exit(status)

    # TODO: Combining classmethod and property works in python 3.9 and later
    # only: https://docs.python.org/3.9/library/functions.html#classmethod
    # Order matters.  In python 3.9, it would be:
    #
    # @classmethod
    # @property
    #
    # Because we're not requiring python 3.9 yet, we have to leave this as a
    # classmethod only:
    @classmethod
    def name(cls):
        """
        Provide the name of the script.  By default, this is the name of the
        module with ``pipnick_`` prepended.
        """
        return f"pipnick_{cls.__module__.split('.')[-1].split('_')[0]}"

    @classmethod
    def parse_args(cls, options=None):
        """
        Parse the command-line arguments.
        """
        parser = cls.get_parser()
        return parser.parse_args() if options is None else parser.parse_args(options)

    # Base classes should override this
    @staticmethod
    def main(args):
        """
        Execute the script.
        """
        pass

    # Base classes should override this.  Ideally they should use this
    # base-class method to intantiate the ArgumentParser object and then fill in
    # the relevant parser arguments
    @classmethod
    def get_parser(
        cls, description=None, width=None, formatter=argparse.ArgumentDefaultsHelpFormatter,
        include_log_options=True
    ):
        """
        Construct the command-line argument parser.

        Parameters
        ----------
        description : str, optional
            A short description of the purpose of the script.
        width : int, optional
            Restrict the width of the formatted help output to be no longer
            than this number of characters, if possible given the help
            formatter.  If None, the width is the same as the terminal width.
        formatter : :class:`argparse.HelpFormatter`, optional
            Class used to format the help output.
        include_log_options : bool, optional
            Include the options controlling logging.

        Returns
        -------
        :class:`argparse.ArgumentParser`
            Command-line interpreter.
        """
        parser = argparse.ArgumentParser(
            description=description, formatter_class=partial(formatter, width=width)
        )
        if not include_log_options:
            return parser
        parser.add_argument(
            '-v', '--verbosity', type=int, default=1,
            help='Verbosity level, which must be 0, 1, or 2.  Level 0 includes warning and '
                 'error messages, level 1 adds informational messages, and level 2 adds '
                 'debugging messages and the calling sequence.'
        )
        parser.add_argument(
            '--log_file', type=str, default=None,
            help='Path for the log file, which must not exist.  By default, the log is written '
                 'to a new file named for the script and the time (UT) in the log directory, '
                 'e.g., <rdxdir>/log.  Logs are never overwritten or appended to.'
        )
        parser.add_argument(
            '--log_level', type=int, default=None,
            help='Verbosity level for the log file.  If a log file is written and this is not '
                 'set, it matches the console verbosity.'
        )
        return parser

    @classmethod
    def init_log(cls, args, log_dir=None):
        """
        Initialize the logger from the parsed command-line arguments, and log
        the pipnick version and the log file.

        The log file is ``args.log_file`` if it is set; otherwise, it is a new
        file in ``log_dir``; see :meth:`default_log_file`.  Logs are never
        overwritten or appended to.

        Parameters
        ----------
        args : :class:`argparse.Namespace`
            The parsed arguments, carrying ``verbosity``, ``log_file`` and
            ``log_level``.
        log_dir : str, :class:`pathlib.Path`, optional
            The directory for the default log file, usually ``<rdxdir>/log``.
            If None and ``args.log_file`` is not set, no log file is written.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickPathError`
            Raised if ``args.log_file`` already exists.
        """
        level = log.convert_verbosity_to_logging_level(args.verbosity)
        log_file_level = (
            None if args.log_level is None
            else log.convert_verbosity_to_logging_level(args.log_level)
        )
        if args.log_file is not None:
            _log_file = Path(args.log_file).absolute()
            if _log_file.exists():
                raise PipnickPathError(
                    f'The log file {_log_file} already exists; logs are never overwritten or '
                    'appended to.  Provide a new file name.'
                )
        elif log_dir is not None:
            _log_file = cls.default_log_file(log_dir)
        else:
            _log_file = None
        log.init(level=level, log_file=_log_file, log_file_level=log_file_level)
        log.info(f'pipnick {version}')
        if _log_file is not None:
            log.info(f'Log file: {_log_file}')

    @staticmethod
    def log_dir(root):
        """
        Return the default log directory.

        Parameters
        ----------
        root : str, :class:`pathlib.Path`
            The directory of the outputs, usually the reduction directory.  If
            None (e.g., it cannot be determined), the current directory is
            used.

        Returns
        -------
        :class:`pathlib.Path`
            ``<root>/log``.
        """
        return (Path.cwd() if root is None else Path(root).absolute()) / 'log'

    @classmethod
    def default_log_file(cls, log_dir):
        """
        Return the default log file: a new file named for the script and the
        current time (UT, to the second), e.g.,
        ``pipnick_reduce_20261028T153614.log``.

        If the file exists (e.g., for two runs started in the same second),
        a counter is appended (``..._2.log``).

        Parameters
        ----------
        log_dir : str, :class:`pathlib.Path`
            The directory for the log file.

        Returns
        -------
        :class:`pathlib.Path`
            The log file, which does not exist.
        """
        timestamp = datetime.datetime.now(datetime.UTC).strftime('%Y%m%dT%H%M%S')
        base = Path(log_dir).absolute() / f'{cls.name()}_{timestamp}.log'
        path = base
        n = 2
        while path.exists():
            path = base.with_name(f'{base.stem}_{n}.log')
            n += 1
        return path
