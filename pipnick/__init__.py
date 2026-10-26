"""
pipnick package initialization.

Provides the package-level globals imported by the submodules.
"""

import logging

from .pkg.version import version
from .pkg.logger import get_logger
# Import all the exceptions so that they can be directly imported (e.g., ``from
# pipnick import PipnickError``).  The names exported are those in
# pipnick.pkg.exceptions.__all__.
from .pkg.exceptions import *  # noqa: F401, F403

# Set version
__version__ = version

# Start the log
log = get_logger(level=logging.INFO)
