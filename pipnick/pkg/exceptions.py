"""
Provides pipnick-specific exceptions.

.. note::

    Adapted from ``pypeit/pkg/exceptions.py`` in `PypeIt
    <https://github.com/pypeit/PypeIt>`__ (BSD 3-Clause).
"""

__all__ = [
    'PipnickError',
    'PipnickCodingError',
    'PipnickDataModelError',
    'PipnickPathError',
    'PipnickParameterError',
]


class PipnickError(Exception):
    """Base class for all pipnick-specific exceptions."""
    pass


class PipnickCodingError(PipnickError):
    """
    Raised for a fault in how pipnick itself is written, not in how it is used.

    This marks a condition that no user input can produce: an invariant the code
    is supposed to maintain, or a branch that should be unreachable.  Raising
    this type rather than :class:`PipnickParameterError` (or similar) keeps the
    two audiences apart -- a user seeing this has found a bug to report, not a
    mistake to correct.
    """
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
