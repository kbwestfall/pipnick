"""
Interfaces to the remote services used for astrometric calibration: the
astrometry.net solver and the Gaia catalog.

These are the only objects in pipnick that contact either service.  Every
method of :class:`AstrometryNetSolver` that takes data submits a job to
astrometry.net (or monitors one), and :meth:`GaiaCatalog.sources` queries the
Gaia archive unless the query is cached.  The pipeline accepts substitutes
with the same methods, which is how it is tested without contacting either
service.

Coordinate frames
-----------------

All astrometric solutions in pipnick are in the ICRS.  Gaia positions are in
the ICRS.  astrometry.net solutions are *labeled* as FK5 (J2000), but are
relabeled as ICRS by :func:`icrs_header`, which explains why this is correct.
astrometry.net provides no setting to request a frame.
"""

import os
from pathlib import Path

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astroquery.astrometry_net import conf as astrometry_net_conf
from astroquery.astrometry_net.core import AstrometryNetClass
from astroquery.exceptions import TimeoutError

from pipnick import log
from pipnick.pkg.exceptions import PipnickError


class AstrometryNetSolver:
    """
    Submit jobs to the astrometry.net web service (https://nova.astrometry.net).

    Parameters
    ----------
    api_key : str, optional
        The astrometry.net API key.  If None, the key is taken from the
        ``ASTROMETRY_NET_API_KEY`` environment variable or, if that is not
        set, from the ``astroquery`` configuration (``api_key`` in the
        ``[astrometry_net]`` section of ``astroquery.cfg``).  The key is never
        logged.
    solve_timeout : int, optional
        Time in seconds to wait for a job to finish before returning it as
        pending; see :meth:`monitor`.

    Notes
    -----
    The solutions are returned as astrometry.net provides them, labeled as
    FK5 (J2000); the pipeline relabels them as ICRS with
    :func:`icrs_header`.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickError`
        Raised if no API key is available.
    """
    def __init__(self, api_key=None, solve_timeout=120):
        key = api_key
        if key is None:
            key = os.environ.get('ASTROMETRY_NET_API_KEY', astrometry_net_conf.api_key)
        if key is None or len(key) == 0:
            raise PipnickError(
                'No astrometry.net API key.  Set the ASTROMETRY_NET_API_KEY environment '
                'variable, or the api_key in the [astrometry_net] section of your astroquery '
                'configuration file.'
            )
        self.solve_timeout = solve_timeout
        self._client = AstrometryNetClass()
        self._client.api_key = key

    @staticmethod
    def _result(header, submission_id):
        """
        Convert a returned solution to the result tuple; see :meth:`solve_sources`.
        """
        if header is None or len(header) == 0:
            return 'failed', None, submission_id
        return 'solved', fits.Header(header), submission_id

    def solve_sources(self, x, y, width, height, hints):
        """
        Submit a list of source positions to astrometry.net and wait for the
        solution.

        Parameters
        ----------
        x, y : :class:`numpy.ndarray`
            The 0-based pixel coordinates of the sources.  They are converted
            to the 1-based (FITS) convention used by astrometry.net.
        width, height : int
            The size of the image in pixels.
        hints : dict
            The settings for the solver; see :func:`astrometry_net_hints`.

        Returns
        -------
        status : str
            ``'solved'``, ``'failed'``, or ``'pending'`` if the job did not
            finish within the timeout.
        header : :class:`astropy.io.fits.Header`, None
            The header with the solution, if solved.
        submission_id : str, None
            The astrometry.net submission ID.
        """
        try:
            header, submission_id = self._client.solve_from_source_list(
                np.asarray(x) + 1, np.asarray(y) + 1, width, height,
                solve_timeout=self.solve_timeout, verbose=False, return_submission_id=True,
                **hints
            )
        except TimeoutError as e:
            return 'pending', None, e.args[1]
        return self._result(header, submission_id)

    def solve_image(self, path, hints):
        """
        Upload an image to astrometry.net, which detects the sources, and wait
        for the solution.

        Parameters
        ----------
        path : :class:`pathlib.Path`
            The FITS file to upload.
        hints : dict
            The settings for the solver; see :func:`astrometry_net_hints`.

        Returns
        -------
        tuple
            The status, header, and submission ID; see :meth:`solve_sources`.
        """
        try:
            header, submission_id = self._client.solve_from_image(
                str(path), force_image_upload=True, solve_timeout=self.solve_timeout,
                verbose=False, return_submission_id=True, **hints
            )
        except TimeoutError as e:
            return 'pending', None, e.args[1]
        return self._result(header, submission_id)

    def monitor(self, submission_id):
        """
        Wait for a pending job to finish.

        Parameters
        ----------
        submission_id : str
            The astrometry.net submission ID.

        Returns
        -------
        tuple
            The status, header, and submission ID; see :meth:`solve_sources`.
        """
        try:
            header, _ = self._client.monitor_submission(
                submission_id, solve_timeout=self.solve_timeout, verbose=False,
                return_submission_id=True
            )
        except TimeoutError:
            return 'pending', None, submission_id
        return self._result(header, submission_id)


def icrs_header(header):
    """
    Relabel an astrometry.net solution as being in the ICRS.

    astrometry.net writes ``EQUINOX = 2000.0`` and no ``RADESYS`` keyword to
    its solutions (see ``util/sip_qfits.c`` in the astrometry.net source), and
    its web API has no setting to request another frame.  By the FITS WCS
    standard, a header with ``EQUINOX`` of 1984 or later and no ``RADESYS``
    is in the FK5 frame, so astropy reads (and writes) these solutions as
    FK5 J2000.

    However, astrometry.net fits its solutions to index files built from
    catalogs whose positions are on the ICRS (Tycho-2, 2MASS, and Gaia), so
    the solutions are effectively in the ICRS; the label is a convention of
    the output, not a property of the solution.  FK5 J2000 and the ICRS
    differ by at most a few tens of milliarcseconds, far below the accuracy
    needed for Nickel data, but treating the solutions as FK5 would mix frames
    in the pipeline: e.g., :class:`astropy.coordinates.SkyCoord` refuses to
    compute offsets between coordinates in different frames.  Relabeling keeps
    every solution in the ICRS, the frame of the Gaia catalog and of the
    telescope pointing used by pipnick.

    Parameters
    ----------
    header : :class:`astropy.io.fits.Header`
        A header with an astrometry.net solution, as returned by the service
        or as cached by pipnick.

    Returns
    -------
    :class:`astropy.io.fits.Header`
        A copy of the header, with ``RADESYS = 'ICRS'`` and without the
        ``EQUINOX`` and (deprecated) ``RADECSYS`` keywords.
    """
    hdr = fits.Header(header).copy()
    for key in ['EQUINOX', 'RADECSYS']:
        hdr.remove(key, ignore_missing=True, remove_all=True)
    hdr['RADESYS'] = ('ICRS', 'astrometry.net solution relabeled as ICRS')
    return hdr


def astrometry_net_hints(center, radius, scale, scale_tolerance=0.1):
    """
    Construct the astrometry.net settings that limit the search for a
    solution.

    Parameters
    ----------
    center : :class:`astropy.coordinates.SkyCoord`
        The expected coordinates of the image center.
    radius : float
        The maximum distance in degrees of the true image center from
        ``center``.
    scale : float
        The expected pixel scale in arcsec per pixel.
    scale_tolerance : float, optional
        The fractional tolerance for the pixel scale.

    Returns
    -------
    dict
        The settings.  The solution is not publicly visible.
    """
    return {
        'center_ra': float(center.icrs.ra.deg), 'center_dec': float(center.icrs.dec.deg),
        'radius': float(radius), 'scale_units': 'arcsecperpix', 'scale_type': 'ul',
        'scale_lower': scale * (1 - scale_tolerance), 'scale_upper': scale * (1 + scale_tolerance),
        'publicly_visible': 'n'
    }


def propagate_proper_motion(catalog, epoch):
    r"""
    Propagate catalog positions to another epoch using their proper motions.

    Parameters
    ----------
    catalog : :class:`astropy.table.Table`
        The catalog, with columns ``ra`` and ``dec`` (deg), ``pmra`` (the
        proper motion in right ascension, including the :math:`\cos\delta`
        factor) and ``pmdec`` (mas/yr), and ``ref_epoch`` (Julian year).
        Masked proper motions are taken to be 0.
    epoch : float
        The epoch (Julian year) of the propagated positions.

    Returns
    -------
    :class:`astropy.coordinates.SkyCoord`
        The propagated coordinates.
    """
    pmra = np.ma.filled(np.ma.asarray(catalog['pmra'], dtype=float), 0.)
    pmdec = np.ma.filled(np.ma.asarray(catalog['pmdec'], dtype=float), 0.)
    dt = epoch - np.asarray(catalog['ref_epoch'], dtype=float)
    dec = np.asarray(catalog['dec'], dtype=float) + pmdec * dt / 3.6e6
    ra = np.asarray(catalog['ra'], dtype=float) \
        + pmra * dt / 3.6e6 / np.cos(np.radians(np.asarray(catalog['dec'], dtype=float)))
    return SkyCoord(ra, dec, unit='deg')


class GaiaCatalog:
    """
    Provide Gaia DR3 sources near a position, caching every query.

    Parameters
    ----------
    cache_dir : str, :class:`pathlib.Path`
        Directory for the cached queries, written as ECSV files.
    mag_limit : float, optional
        The faintest G-band magnitude to return.
    max_sources : int, optional
        The maximum number of (brightest) sources returned by each query.
    pad : float, optional
        Queries are made for a radius this factor larger than requested, so
        that nearby images can reuse them.
    """
    def __init__(self, cache_dir, mag_limit=19., max_sources=1000, pad=1.5):
        self.cache_dir = Path(cache_dir)
        self.mag_limit = mag_limit
        self.max_sources = max_sources
        self.pad = pad

    def _find_cache(self, center, radius):
        """
        Find a cached query that covers a cone, or return None.
        """
        if not self.cache_dir.is_dir():
            return None
        for f in sorted(self.cache_dir.glob('gaia_*.ecsv')):
            tbl = Table.read(f, format='ascii.ecsv')
            cached = SkyCoord(tbl.meta['center_ra'], tbl.meta['center_dec'], unit='deg')
            if (
                tbl.meta['mag_limit'] >= self.mag_limit
                and cached.separation(center).deg + radius <= tbl.meta['radius']
            ):
                log.info(f'Using the cached Gaia query in {f.name}.')
                return tbl
        return None

    def _query(self, center, radius):
        """
        Query the Gaia archive.  This is the only method that contacts it.
        """
        # NOTE: astroquery.gaia contacts the Gaia archive when it is imported, so
        # it is only imported when a query is actually made.
        from astroquery.gaia import Gaia
        query = (
            f'SELECT TOP {self.max_sources} ra, dec, pmra, pmdec, ref_epoch, phot_g_mean_mag '
            'FROM gaiadr3.gaia_source WHERE 1 = CONTAINS(POINT(ra, dec), '
            f'CIRCLE({center.icrs.ra.deg}, {center.icrs.dec.deg}, {radius})) '
            f'AND phot_g_mean_mag < {self.mag_limit} ORDER BY phot_g_mean_mag'
        )
        log.info(
            f'Querying Gaia DR3 within {radius*60:.1f} arcmin of {center.to_string("hmsdms")}.'
        )
        return Table(Gaia.launch_job(query).get_results())

    def sources(self, center, radius, epoch):
        """
        Return the Gaia sources within a cone.

        Parameters
        ----------
        center : :class:`astropy.coordinates.SkyCoord`
            The center of the cone.
        radius : float
            The radius of the cone in degrees.
        epoch : float
            The epoch (Julian year) of the returned positions.

        Returns
        -------
        coords : :class:`astropy.coordinates.SkyCoord`
            The coordinates of the sources, propagated to ``epoch``; see
            :func:`propagate_proper_motion`.
        mag : :class:`numpy.ndarray`
            The G-band magnitudes.
        """
        tbl = self._find_cache(center, radius)
        if tbl is None:
            qradius = radius * self.pad
            tbl = self._query(center, qradius)
            tbl.meta = {
                'center_ra': float(center.icrs.ra.deg), 'center_dec': float(center.icrs.dec.deg),
                'radius': float(qradius), 'mag_limit': float(self.mag_limit)
            }
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            ofile = self.cache_dir / (
                f'gaia_{center.icrs.ra.deg:08.4f}_{center.icrs.dec.deg:+08.4f}_'
                f'{qradius*60:.1f}.ecsv'
            )
            # Only the needed columns, with simple names
            tbl = Table(
                {
                    'ra': np.asarray(tbl['ra'], dtype=float),
                    'dec': np.asarray(tbl['dec'], dtype=float),
                    'pmra': np.ma.asarray(tbl['pmra'], dtype=float),
                    'pmdec': np.ma.asarray(tbl['pmdec'], dtype=float),
                    'ref_epoch': np.asarray(tbl['ref_epoch'], dtype=float),
                    'mag': np.asarray(tbl['phot_g_mean_mag'], dtype=float)
                },
                meta=tbl.meta
            )
            tbl.write(ofile, format='ascii.ecsv', overwrite=True)
            log.info(f'Cached the Gaia query in {ofile.name}.')
        coords = propagate_proper_motion(tbl, epoch)
        inside = coords.separation(center).deg <= radius
        return coords[inside], np.asarray(tbl['mag'], dtype=float)[inside]
