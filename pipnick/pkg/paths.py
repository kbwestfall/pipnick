"""
The names of the products that pipnick writes in a reduction directory.

The reduced images are written to the reduction directory itself (see
:func:`~pipnick.pipelines.reduction.reduce_all`).  The products of each later
step are written to a subdirectory named for the step (see :func:`step_dir`),
and the products derived from one reduced image are named for the image.
Steps that read the products of another step use these functions to find
them.
"""

from pathlib import Path

from pipnick.pkg.exceptions import PipnickParameterError


def step_dir(rdxdir, step):
    """
    Return the directory for the products of a step.

    Parameters
    ----------
    rdxdir : str, :class:`pathlib.Path`
        The reduction directory.
    step : str
        The step: ``'astrometry'``, ``'photometry'``, ``'calibrate'`` or
        ``'mosaic'``.  The directory has the same name.

    Returns
    -------
    :class:`pathlib.Path`
        The directory.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickParameterError`
        Raised if the step is not one of the above.
    """
    if step not in ('astrometry', 'photometry', 'calibrate', 'mosaic'):
        raise PipnickParameterError(f'Unknown pipnick step: {step}.')
    return Path(rdxdir) / step


def photometry_table(image):
    """
    Return the source table of a reduced image; see
    :func:`~pipnick.pipelines.photometry.photometry_all`.

    Parameters
    ----------
    image : str, :class:`pathlib.Path`
        The reduced image, in the reduction directory.

    Returns
    -------
    :class:`pathlib.Path`
        The table, ``photometry/<image>_phot.ecsv`` in the reduction directory.
    """
    _image = Path(image)
    return step_dir(_image.parent, 'photometry') / f'{_image.stem}_phot.ecsv'


def points_file(image):
    """
    Return the file with the astrometric solution of a reduced image provided
    by the user; see :func:`~pipnick.pipelines.astrometry.apply_wcs`.

    Parameters
    ----------
    image : str, :class:`pathlib.Path`
        The reduced image, in the reduction directory.

    Returns
    -------
    :class:`pathlib.Path`
        The file, ``astrometry/<image>_points.fits`` in the reduction
        directory.
    """
    _image = Path(image)
    return step_dir(_image.parent, 'astrometry') / f'{_image.stem}_points.fits'


def calibrated_table(image):
    """
    Return the calibrated source table of a reduced image; see
    :func:`~pipnick.pipelines.calibrate.apply_to_catalogs`.

    Parameters
    ----------
    image : str, :class:`pathlib.Path`
        The reduced image, in the reduction directory.

    Returns
    -------
    :class:`pathlib.Path`
        The table, ``calibrate/<image>_cal.ecsv`` in the reduction directory.
    """
    _image = Path(image)
    return step_dir(_image.parent, 'calibrate') / f'{_image.stem}_cal.ecsv'


def calibration_file(rdxdir, tag, night):
    """
    Return the photometric calibration file of a night; see
    :func:`~pipnick.pipelines.calibrate.calibrate_all`.

    Parameters
    ----------
    rdxdir : str, :class:`pathlib.Path`
        The reduction directory.
    tag : str
        The setup tag; see :func:`~pipnick.metadata.setup_tag`.
    night : str
        The night (UT date, ``YYYY-MM-DD``).

    Returns
    -------
    :class:`pathlib.Path`
        The file, ``calibrate/<tag>_<night>.toml`` in the reduction directory.
    """
    return step_dir(rdxdir, 'calibrate') / f'{tag}_{night}.toml'


def residuals_file(rdxdir, tag, night):
    """
    Return the table of the standard-star measurements used for the
    photometric calibration of a night; see
    :func:`~pipnick.pipelines.calibrate.calibrate_all`.

    Parameters
    ----------
    rdxdir : str, :class:`pathlib.Path`
        The reduction directory.
    tag : str
        The setup tag; see :func:`~pipnick.metadata.setup_tag`.
    night : str
        The night (UT date, ``YYYY-MM-DD``).

    Returns
    -------
    :class:`pathlib.Path`
        The table, ``calibrate/<tag>_<night>_residuals.ecsv`` in the
        reduction directory.
    """
    return step_dir(rdxdir, 'calibrate') / f'{tag}_{night}_residuals.ecsv'
