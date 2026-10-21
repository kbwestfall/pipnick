"""
pipnick package initialization.

Provides the package-level globals imported by the submodules.
"""

from .pkg.version import version

# Set version
__version__ = version

# Start the log
import logging
from .pkg.logger import get_logger  # noqa: E402
log = get_logger(level=logging.INFO)

# Import all the exceptions so that they can be directly imported (e.g., `from
# pipnick import PipnickError`) in all package imports.  The names exported are
# those in pipnick.pkg.exceptions.__all__.
from .pkg.exceptions import *  # noqa: E402, F401, F403
