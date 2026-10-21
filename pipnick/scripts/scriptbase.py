"""
Implements base classes for use with scripts.

.. note::

    Adapted from ``pypeit/scripts/scriptbase.py`` in `PypeIt
    <https://github.com/pypeit/PypeIt>`__ (BSD 3-Clause); see
    ``licenses/README.rst``.
"""

import argparse
import datetime
import textwrap
from functools import reduce

from .. import log

class SmartFormatter(argparse.HelpFormatter):
    r"""
    Enable a combination of both fixed-format and wrappable lines to be
    formatted for the help statements for command-line arguments used with
    :class:`argparse.ArgumentParser`.

    Borrows from
    https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text

    Help strings that use this formatter *must* begin with "R|".  If not, the
    help string is parsed by the base class.

    When parsed by this formatter, the leading "R|" characters are stripped and
    the lines to be printed are parsed using :meth:`str.splitlines`.  Each resulting
    line is wrapped using :func:`textwrap.wrap`, unless it begins with the characters
    "F|", which forces the line to remain unaltered (except for stripping the
    leading characters).

    For example, if you add an argument like this:

    .. code-block:: python

        parser.add_argument('-t', '--tell_file', type=str,
                            help='R|Configuration file to change default telluric parameters.  '
                                 'Note that the parameters in this file will be overwritten if '
                                 'you set argument in your terminal.  The --tell_file option '
                                 'requires a .tell file with the following format:\n'
                                 '\n'
                                 'F|    [tellfit]\n'
                                 'F|         objmodel = qso\n'
                                 'F|         redshift = 7.6\n'
                                 'F|         bal_wv_min_max = 10825,12060\n'
                                 'OR\n'
                                 'F|    [tellfit]\n'
                                 'F|         objmodel = star\n'
                                 'F|         star_type = A0\n'
                                 'F|         star_mag = 8.\n'
                                 'OR\n'
                                 'F|    [tellfit]\n'
                                 'F|         objmodel = poly\n'
                                 'F|         polyorder = 3\n'
                                 'F|         fit_wv_min_max = 9000.,9500.\n'
                                 '\n')

    The result will be (depending on the width of your console):

    .. code-block:: console

        -t TELL_FILE, --tell_file TELL_FILE
                          Configuration file to change default telluric
                          parameters.  Note that the parameters in this file
                          will be overwritten if you set argument in your
                          terminal.  The --tell_file option requires a .tell
                          file with the following format:

                              [tellfit]
                                   objmodel = qso
                                   redshift = 7.6
                                   bal_wv_min_max = 10825,12060
                          OR
                              [tellfit]
                                   objmodel = star
                                   star_type = A0
                                   star_mag = 8.
                          OR
                              [tellfit]
                                   objmodel = poly
                                   polyorder = 3
                                   fit_wv_min_max = 9000.,9500.
    """
    def _split_lines(self, text, width):
        """
        Split the provided text into width constrained lines.

        See the class description for formatting instructions.
        """
        if text.startswith('R|'):
            lines = text[2:].splitlines()
            for i in range(len(lines)):
                if lines[i].startswith('F|'):
                    lines[i] = [lines[i][2:]]
                elif len(lines[i]) == 0:
                    lines[i] = [' ']
                else:
                    lines[i] = textwrap.wrap(lines[i], width)
            return reduce(list.__add__, lines)
        return super()._split_lines(text, width)


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
    def get_parser(cls, description=None, width=None,
                   formatter=argparse.ArgumentDefaultsHelpFormatter, include_log_options=True):
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
        formatter : argparse.HelpFormatter, optional
            Class used to format the help output.
        include_log_options : bool, optional
            Include the options controlling logging.

        Returns
        -------
        :class:`argparse.ArgumentParser`
            Command-line interpreter.
        """
        parser = argparse.ArgumentParser(description=description,
                                         formatter_class=lambda prog: formatter(prog, width=width))
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


