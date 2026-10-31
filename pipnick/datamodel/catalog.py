"""
The source table of a reduced image.

The photometry writes the sources measured in each reduced image to a table,
``photometry/<image>_phot.ecsv`` (see
:func:`~pipnick.pipelines.photometry.photometry_all`).  :class:`SourceCatalog`
wraps this table; its columns and metadata are accessed as those of the table
(e.g., ``catalog['x']`` and ``catalog.meta``).
"""

import numpy as np

from astropy.table import MaskedColumn, Table

from pipnick.astrometry.fit import to_icrs, wcs_from_header
from pipnick.datamodel.psf import PSFModel


def celestial_coordinates(header, x, y):
    """
    Compute the celestial coordinates of sources, if the image has an
    astrometric solution found by pipnick.

    Parameters
    ----------
    header : :class:`astropy.io.fits.Header`
        The header of the reduced image.
    x, y : :class:`numpy.ndarray`
        The 0-based pixel coordinates of the sources.

    Returns
    -------
    ra, dec : :class:`astropy.table.MaskedColumn`
        The ICRS right ascension and declination in degrees, masked if the
        image has no solution.
    """
    n = len(x)
    if 'PNASTMTH' not in header or n == 0:
        masked = np.ones(n, dtype=bool)
        return (
            MaskedColumn(np.zeros(n), name='ra', unit='deg', mask=masked),
            MaskedColumn(np.zeros(n), name='dec', unit='deg', mask=masked)
        )
    coords = to_icrs(wcs_from_header(header).pixel_to_world(x, y))
    unmasked = np.zeros(n, dtype=bool)
    return (
        MaskedColumn(coords.ra.deg, name='ra', unit='deg', mask=unmasked),
        MaskedColumn(coords.dec.deg, name='dec', unit='deg', mask=unmasked)
    )


class SourceCatalog:
    """
    The source table of a reduced image.

    Parameters
    ----------
    table : :class:`astropy.table.Table`
        The table; see :func:`~pipnick.pipelines.photometry.photometry_all`.
    path : :class:`pathlib.Path`, optional
        The file of the table.
    """
    def __init__(self, table, path=None):
        self.table = table
        self.path = path

    @classmethod
    def read(cls, path):
        """
        Read a source table.

        Parameters
        ----------
        path : :class:`pathlib.Path`
            The file.

        Returns
        -------
        :class:`SourceCatalog`
            The catalog.
        """
        return cls(Table.read(path, format='ascii.ecsv'), path=path)

    def write(self, path=None):
        """
        Write the source table.

        Parameters
        ----------
        path : :class:`pathlib.Path`, optional
            The file; if None, the file it was read from.
        """
        _path = self.path if path is None else path
        _path.parent.mkdir(parents=True, exist_ok=True)
        self.table.write(_path, format='ascii.ecsv', overwrite=True)

    def __len__(self):
        return len(self.table)

    def __getitem__(self, key):
        return self.table[key]

    @property
    def meta(self):
        """The metadata of the table."""
        return self.table.meta

    @property
    def colnames(self):
        """The column names of the table."""
        return self.table.colnames

    @property
    def psf(self):
        """The PSF model of the image, or None if there is none."""
        return None if 'psf' not in self.meta else PSFModel.from_meta(self.meta['psf'])

    @property
    def has_coordinates(self):
        """True if the sources have celestial coordinates; see :meth:`update_coordinates`."""
        return bool(
            self.meta.get('wcs_method', 'none') != 'none'
            and np.any(np.logical_not(np.ma.getmaskarray(np.ma.asarray(self.table['ra']))))
        )

    def update_coordinates(self, header, wcs_hash):
        """
        Update the celestial coordinates of the sources for a new astrometric
        solution; see :func:`celestial_coordinates`.

        Parameters
        ----------
        header : :class:`astropy.io.fits.Header`
            The header of the reduced image.
        wcs_hash : str
            The hash of the solution; see
            :func:`~pipnick.pkg.provenance.solution_hash`.
        """
        self.table['ra'], self.table['dec'] = celestial_coordinates(
            header, np.asarray(self.table['x']), np.asarray(self.table['y'])
        )
        self.meta['wcs_method'] = header.get('PNASTMTH', 'none')
        self.meta['wcs_hash'] = wcs_hash

    def usable(self):
        """
        Select the sources usable for relative photometry: those that are not
        flagged, have an unmasked aperture flux (so no saturated or masked
        pixels), and have celestial coordinates.

        Returns
        -------
        :class:`numpy.ndarray`
            Boolean array selecting the sources.
        """
        return (
            (np.asarray(self.table['flags']) == 0)
            & np.logical_not(np.ma.getmaskarray(np.ma.asarray(self.table['flux_aper'])))
            & np.logical_not(np.ma.getmaskarray(np.ma.asarray(self.table['ra'])))
        )

    def instrumental_mag(self, mag):
        r"""
        Return the instrumental magnitudes of the sources and their errors.

        For ``mag='psf'``, these are ``mag_inst`` and ``mag_inst_err``.  For
        ``mag='aper'``, they are :math:`-2.5\log(f / t)` for the aperture flux
        :math:`f` and the exposure time :math:`t`, masked where the flux is
        masked or not positive.

        Parameters
        ----------
        mag : str
            ``'psf'`` or ``'aper'``.

        Returns
        -------
        mag, mag_err : :class:`numpy.ma.MaskedArray`
            The magnitudes and their errors.
        """
        if mag == 'psf':
            return (
                np.ma.asarray(self.table['mag_inst'], dtype=float),
                np.ma.asarray(self.table['mag_inst_err'], dtype=float)
            )
        flux = np.ma.asarray(self.table['flux_aper'], dtype=float)
        ferr = np.ma.asarray(self.table['flux_aper_err'], dtype=float)
        bad = np.ma.getmaskarray(flux) | (flux.filled(0.) <= 0)
        _flux = np.where(bad, 1., flux.filled(1.))
        value = -2.5 * np.log10(_flux / float(self.meta['exptime']))
        err = 2.5 / np.log(10.) * ferr.filled(0.) / _flux
        return (
            np.ma.MaskedArray(value, mask=bad),
            np.ma.MaskedArray(err, mask=bad | np.ma.getmaskarray(ferr))
        )

    def neighbors(self, index, residual_radius=1.5, max_neighbor_flux=0.01):
        """
        Classify the sources near a source.

        The PSF photometry of bright stars often detects faint, spurious sources
        in the residuals of their fits, very close to the star, which are then fit
        together with it (so its ``group_size`` is larger than 1).  These
        *residual detections* are identified as sources detected in a later
        iteration than the source (``iter_detected``) and within
        ``residual_radius`` FWHM of it; they are counted, but do not exclude the
        source.  Any other source within the photometric aperture of the source
        (``aper_radius_pix`` in the table metadata) with more than
        ``max_neighbor_flux`` of its PSF flux is a *contaminating* neighbor.

        Parameters
        ----------
        index : int
            The row of the source.
        residual_radius : float, optional
            The radius, in FWHM, within which later detections are taken to be in
            the residuals of the source.
        max_neighbor_flux : float, optional
            The fraction of the PSF flux of the source above which a neighbor
            within the aperture contaminates it.

        Returns
        -------
        n_residual, n_contaminating : int
            The numbers of residual detections and contaminating neighbors.
        """
        x = np.asarray(self.table['x'], dtype=float)
        y = np.asarray(self.table['y'], dtype=float)
        dist = np.hypot(x - x[index], y - y[index])
        flux = np.ma.filled(np.ma.asarray(self.table['flux_psf'], dtype=float), 0.)
        iteration = (
            np.asarray(self.table['iter_detected']) if 'iter_detected' in self.table.colnames
            else np.ones(len(self.table), dtype=int)
        )
        other = np.arange(len(self.table)) != index
        residual = (
            other & (iteration > iteration[index])
            & (dist < residual_radius * self.meta['psf']['fwhm_pix'])
        )
        contaminating = (
            other & np.logical_not(residual) & (dist < self.meta['aper_radius_pix'])
            & (flux > max_neighbor_flux * flux[index])
        )
        return int(np.sum(residual)), int(np.sum(contaminating))
