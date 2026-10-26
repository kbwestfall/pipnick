"""
Provides pipnick-specific exceptions.

.. note::

    Adapted from ``pypeit/pkg/exceptions.py`` in `PypeIt
    <https://github.com/pypeit/PypeIt>`__ (BSD 3-Clause); see
    ``licenses/README.rst``.
"""

__all__ = ['PipnickError', 'PipnickDataModelError', 'PipnickPathError', 'PipnickParameterError']


class PipnickError(Exception):
    """Base class for all pipnick-specific exceptions."""
    pass


class PipnickDataModelError(PipnickError):
    """
    Raised when data (e.g., a metadata table or FITS header) does not conform
    to its expected datamodel.
    """
    pass


class PipnickPathError(PipnickError):
    """Raised for missing or invalid filesystem paths."""
    pass


class PipnickParameterError(PipnickError):
    """Raised for invalid parameter values."""
    pass
