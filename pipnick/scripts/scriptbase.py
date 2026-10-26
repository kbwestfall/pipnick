"""
Implements base classes for use with scripts.

.. note::

    Adapted from ``pypeit/scripts/scriptbase.py`` in `PypeIt
    <https://github.com/pypeit/PypeIt>`__ (BSD 3-Clause); see
    ``licenses/README.rst``.
"""

import argparse
import datetime
from functools import partial

from .. import log


class ScriptBase:
    """
    Provides a base class for all scripts.
    """
    @classmethod
    def entry_point(cls):
        """
        Defines the main script entry point.
        """
        cls.main(cls.parse_args())

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
            help='Name for the log file.  If set to "default", a name is generated from the '
                 'script name and the current time.  If not set, no log file is written.'
        )
        parser.add_argument(
            '--log_level', type=int, default=None,
            help='Verbosity level for the log file.  If a log file is written and this is not '
                 'set, it matches the console verbosity.'
        )
        return parser

    @classmethod
    def init_log(cls, args):
        """
        Initialize the logger from the parsed command-line arguments.

        Parameters
        ----------
        args : :class:`argparse.Namespace`
            The parsed arguments, carrying ``verbosity``, ``log_file`` and
            ``log_level``.
        """
        level = log.convert_verbosity_to_logging_level(args.verbosity)
        log_file_level = (
            None if args.log_level is None
            else log.convert_verbosity_to_logging_level(args.log_level)
        )
        if args.log_file == 'default':
            _log_file = cls.default_log_file()
        elif args.log_file in ['None', None]:
            _log_file = None
        else:
            _log_file = args.log_file
        log.init(level=level, log_file=_log_file, log_file_level=log_file_level)

    @classmethod
    def default_log_file(cls):
        """
        Return the default name for the log file.

        Returns
        -------
        str
            The script name and a UT timestamp, to the minute.
        """
        timestamp = datetime.datetime.now(datetime.UTC).strftime('%Y%m%d-%H%M')
        return f'{cls.name()}_{timestamp}.log'


