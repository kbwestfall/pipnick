"""
A reduced science image and the products derived from it.
"""

import numpy as np

from astropy.io import fits
from astropy.time import Time
from astropy.wcs.utils import proj_plane_pixel_area

from pipnick.astrometry.fit import wcs_from_header
from pipnick.datamodel.catalog import SourceCatalog
from pipnick.image.io import read_reduced
from pipnick.pkg.exceptions import PipnickDataModelError
from pipnick.pkg.paths import photometry_table
from pipnick.pkg.provenance import image_hash, names_hash, solution_hash, values_hash


class ReducedImage:
    """
    A reduced science image of a data set.

    Reduced images are usually constructed by
    :meth:`~pipnick.datamodel.dataset.DataSet.science_images`.  The header is
    read when first needed; use :meth:`reload` after the file changes.

    Parameters
    ----------
    product : :class:`~pipnick.datamodel.dataset.Product`
        The product of the image.
    dataset : :class:`~pipnick.datamodel.dataset.DataSet`
        The data set of the image.
    """
    def __init__(self, product, dataset):
        self.product = product
        self.dataset = dataset
        self.path = product.file
        self._header = None
        self._values = None
        self._phot = None

    @property
    def name(self):
        """The file name of the image."""
        return self.path.name

    @property
    def exists(self):
        """True if the image has been reduced."""
        return self.path.is_file()

    @property
    def header(self):
        """The primary header of the image."""
        if self._header is None:
            self._header = fits.getheader(self.path)
        return self._header

    def reload(self):
        """Forget the header and source table read so far, e.g., after the file changes."""
        self._header = None
        self._phot = None

    def read(self):
        """
        Read the image; see :func:`~pipnick.image.io.read_reduced`.

        Returns
        -------
        data, mask, error, header : tuple
            The image, its mask and error, and its header.
        """
        data, mask, error, header = read_reduced(self.path)
        self._header = header
        return data, mask, error, header

    def reduction_status(self):
        """
        Determine if the image is up to date with its raw frames, metadata
        table values, and calibrations, as the reduction does (see
        :func:`~pipnick.pipelines.reduction.reduce_all`).

        A calibration that has not been reduced is expected to be absent
        (``'none'``) from the provenance of the image.

        Returns
        -------
        status : str
            ``'ok'``, ``'stale'`` (the image must be reduced again), or
            ``'missing'`` (the image has not been reduced).
        reason : str
            The reason for a stale or missing image.
        """
        if not self.exists:
            return 'missing', 'not reduced'
        hdr = self.header
        if hdr.get('PNRAWHSH') != names_hash(self.dataset.raw_files(self.product)):
            return 'stale', 'the raw frames changed'
        if hdr.get('PNTABH') != values_hash(self.values):
            return 'stale', 'the metadata table values changed'
        for name, calib in self.dataset.calibrations(self.product).items():
            expected = (
                'none' if calib is None or not calib.file.is_file()
                else names_hash(self.dataset.raw_files(calib))
            )
            if hdr.get(f'PN{name.upper()}H') != expected:
                return 'stale', f'the {name} changed'
        return 'ok', '-'

    def check_setup(self):
        """
        Check that the image was reduced from its data set.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
            Raised if the setup of the image (``PNSETUP``) is not that of its
            data set.
        """
        if self.header.get('PNSETUP') != self.dataset.tag:
            raise PipnickDataModelError(f'{self.name} is not from setup {self.dataset.tag}.')

    @property
    def values(self):
        """
        The metadata table values of the image; see
        :meth:`~pipnick.datamodel.dataset.DataSet.table_values`.
        """
        if self._values is None:
            self._values = self.dataset.table_values(self.product)
        return self._values

    @property
    def filter(self):
        """The filter, from the metadata table (None if unknown)."""
        return self.values['filter']

    @property
    def exptime(self):
        """The exposure time (s), from the metadata table (None if unknown)."""
        return self.values['exptime']

    @property
    def airmass(self):
        """The airmass, from the metadata table (None if unknown)."""
        return self.values['airmass']

    @property
    def mjd(self):
        """The MJD of the first frame (None if unknown)."""
        return self.dataset.mjd(self.product)

    @property
    def epoch(self):
        """The epoch of the first frame (Julian year; None if unknown)."""
        mjd = self.mjd
        return None if mjd is None else Time(mjd, format='mjd').jyear

    @property
    def night(self):
        """The UT date of the first frame (``YYYY-MM-DD``; None if unknown)."""
        return self.dataset.night(self.product)

    @property
    def pointing(self):
        """The telescope pointing, from the metadata table (None if unknown)."""
        return self.dataset.pointing(self.product)

    @property
    def nominal_scale(self):
        """
        The nominal pixel scale (arcsec per pixel), from the camera and the
        binning; see :meth:`~pipnick.datamodel.dataset.DataSet.pixel_scale`.
        """
        return self.dataset.pixel_scale(self.product)

    @property
    def shape(self):
        """The shape of the image (rows, columns)."""
        return (int(self.header['NAXIS2']), int(self.header['NAXIS1']))

    @property
    def has_solution(self):
        """True if the image has an astrometric solution found by pipnick."""
        return 'PNASTMTH' in self.header

    @property
    def wcs(self):
        """The WCS of the image; see :func:`~pipnick.astrometry.fit.wcs_from_header`."""
        return wcs_from_header(self.header)

    @property
    def pixel_scale(self):
        """The pixel scale (arcsec per pixel) of the WCS of the image."""
        return float(np.sqrt(proj_plane_pixel_area(self.wcs)) * 3600)

    @property
    def image_hash(self):
        """The hash of the image content; see :func:`~pipnick.pkg.provenance.image_hash`."""
        return image_hash(self.header)

    @property
    def solution_hash(self):
        """
        The hash of the astrometric solution of the image; see
        :func:`~pipnick.pkg.provenance.solution_hash`.
        """
        return solution_hash(self.header)

    @property
    def zero_point(self):
        """The photometric zero point (``PNCALZP``), or None if not calibrated."""
        return self.header.get('PNCALZP')

    @property
    def phot_file(self):
        """The source table of the image; see :func:`~pipnick.pkg.paths.photometry_table`."""
        return photometry_table(self.path)

    def photometry(self, require=True):
        """
        Return the source table of the image, after checking that it is up to
        date.

        Parameters
        ----------
        require : bool, optional
            Raise an exception if the image has no source table.  Otherwise,
            None is returned.

        Returns
        -------
        :class:`~pipnick.datamodel.catalog.SourceCatalog`
            The source table; see
            :func:`~pipnick.pipelines.photometry.photometry_all`.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
            Raised if the table is required but missing, or is out of date:
            the image or its astrometric solution has changed since it was
            written.
        """
        if self._phot is not None:
            return self._phot
        pfile = self.phot_file
        if not pfile.is_file():
            if require:
                raise PipnickDataModelError(
                    f'{self.name} has no source table; run pipnick_photometry.'
                )
            return None
        phot = SourceCatalog.read(pfile)
        if phot.meta.get('image_hash') != self.image_hash:
            raise PipnickDataModelError(
                f'The source table {pfile.name} is out of date: {self.name} has changed since '
                'it was measured.  Rerun pipnick_photometry.'
            )
        if phot.meta.get('wcs_hash') != self.solution_hash:
            raise PipnickDataModelError(
                f'The celestial coordinates in {pfile.name} are out of date: the astrometric '
                f'solution of {self.name} has changed.  Rerun pipnick_photometry to update '
                'them.'
            )
        self._phot = phot
        return phot
