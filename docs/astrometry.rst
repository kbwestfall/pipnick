Astrometric Solution
====================

``pipnick_astrometry`` finds the astrometric solution (WCS) of each reduced
science image and writes it to the image header.  It uses the
`astrometry.net <https://nova.astrometry.net>`__ web service, but submits as
few images as possible: most images are solved locally, by matching their
sources to an image that is already solved or to the Gaia catalog.

Overview
--------

The images are those produced by ``pipnick_reduce`` from a per-setup metadata
file (see ``pipnick_metadata --split``).  For each run:

1. **Existing solutions are kept.**  Solutions that were lost because an image
   was rewritten (e.g., by rerunning ``pipnick_reduce``) are restored from the
   cache of earlier astrometry.net results, or from solutions you provided
   (see `Providing a solution`_), without submitting anything.

2. **Images are grouped by pointing.**  Images whose telescope pointings are
   within ``--group_radius`` (default 1 arcmin) form a group.  One image per
   group, the reference, is solved, and the others are solved by matching
   their detected sources to those of the reference (*bootstrapping*).

3. **Groups without a solved image** are solved from an overlapping solved
   image, or by matching their sources to the Gaia DR3 catalog.  The match
   starts from the solution predicted from a solved image in the same setup
   and the telescope pointing offset measured from the solved images.

4. **Images that still cannot be solved are matched to Gaia with a wide
   search**, which allows for an incorrect telescope pointing in the
   metadata table (see `Incorrect telescope pointings`_).

5. **Only images that cannot be solved this way are submitted** to
   astrometry.net, starting with the image with the most detected sources.
   Typically, one submission is enough for a night.

A solution found by bootstrapping or from Gaia is accepted only if at least
``--min_matches`` sources (default 4) are matched and the RMS residual is at
most ``--max_rms`` (default 0.3 arcsec).  Images are bootstrapped from the
solved image that overlaps them best, preferring images that were not
themselves bootstrapped.

By default, a submission sends only a list of sources detected locally (the
positions of up to 100 of the brightest sources), not the image.  The
solutions are not publicly visible on astrometry.net.

The astrometry.net API key
--------------------------

Register at https://nova.astrometry.net and copy your API key from your
profile.  ``pipnick`` reads it from the ``ASTROMETRY_NET_API_KEY`` environment
variable (e.g., set in your shell startup file) or, if that is not set, from
the ``api_key`` entry in the ``[astrometry_net]`` section of your
``astroquery`` configuration file.  The key is never a command-line argument
and is never logged.

Using the command-line script
-----------------------------

.. code-block:: console

   pipnick_astrometry -t n1_2x2_slow_full.ecsv

The options are:

- ``-t``, ``--raw_table``: The per-setup metadata file used to reduce the
  images (required).
- ``-r``, ``--rdxdir``: The reduction directory.  If not provided, the
  directory saved in the metadata file is used.
- ``--dry_run``: Detect the sources, group the images, and report what would
  be submitted, without solving or submitting anything.  The source lists
  that would be submitted are written to ``astrometry/dry_run``.
- ``--max_submissions``: The maximum number of submissions to astrometry.net.
- ``--resolve``: Ignore the existing and cached solutions and solve every image
  again.  This can submit images that were solved before.
- ``--upload_image``: Upload the images (with masked pixels set to 0) instead
  of source lists, so that astrometry.net detects the sources.
- ``--group_radius``, ``--min_matches``, ``--max_rms``: See `Overview`_.
- ``--wide_search``, ``--wide_min_matches``: The maximum pointing error in
  arcmin allowed by the wide search (default 10; 0 disables it), and the
  minimum number of matched sources it requires (default 10).  See
  `Incorrect telescope pointings`_.

A good first step for a new night is a dry run, followed by a run limited to
one submission:

.. code-block:: console

   pipnick_astrometry -t n1_2x2_slow_full.ecsv --dry_run
   pipnick_astrometry -t n1_2x2_slow_full.ecsv --max_submissions 1

Rerunning the script is safe: solved images are not changed, and nothing is
submitted unless an image cannot be solved otherwise.

Outputs
-------

The solution is written to the primary header of each reduced image,
replacing any existing WCS, with these keywords recording how it was found:

- ``PNASTMTH``: The method: ``astrometry.net``, ``bootstrap`` (from another
  image), ``gaia``, or ``points`` (provided by you).
- ``PNASTWH``: A hash of the solution.
- ``PNASTNPT`` and ``PNASTRMS``: The number of matched sources and the RMS
  residual in arcsec (for bootstrapped, Gaia, and provided solutions).
- ``PNASTREF`` and ``PNASTRH``: The reference image (or ``Gaia DR3``) and the
  hash of its solution.  If the solution of a reference image changes, the
  images bootstrapped from it are solved again on the next run.
- ``PNASTIN`` and ``PNASTSUB``: The hash of the submitted input and the
  astrometry.net submission ID.
- ``PNASTPOE``: For images solved by the wide search, the distance in arcmin
  between the image center and the pointing in the metadata table.

The ``astrometry`` subdirectory of the reduction directory holds the cached
astrometry.net results (``anet_<hash>.fits``, including failed and pending
submissions, which are not resubmitted), the cached Gaia queries
(``gaia_*.ecsv``), and the solutions you provided (``<image>_points.fits``).
Do not delete it unless you want to submit the images again.

Incorrect telescope pointings
-----------------------------

The pointing (RA and Dec) of each image is taken from the metadata table,
which takes it from the raw header.  Occasionally, the raw header is wrong;
e.g., the telescope coordinates may not have been updated after a slew, so
that a series of images carry the pointing of an earlier field.  The images
are then not where ``pipnick_astrometry`` expects them, and the usual
matching fails.

There are two ways to solve such images, neither of which requires a
submission to astrometry.net:

- **Automatically, with the wide search.**  An image that cannot be solved
  otherwise is matched to Gaia allowing for a pointing error of up to
  ``--wide_search`` arcmin (default 10), as long as another image in the
  same setup is solved.  Because chance matches are more likely over a wide
  area, at least ``--wide_min_matches`` sources (default 10) must match.  The
  pointing error is recorded in the header (``PNASTPOE``) and reported with a
  warning that gives the true image center.  Other images with the same
  (stale) pointing then use the measured error, so they are solved by the
  usual matching, even if they have too few sources for the wide search.

- **By correcting the metadata table.**  Edit the ``ra`` and ``dec`` of the
  affected frames in the metadata file (e.g., using the image centers
  reported by the wide search, or the coordinates in the observing log), and
  rerun the reduction and the astrometry:

  .. code-block:: console

     pipnick_reduce -t n1_2x2_slow_full.ecsv
     pipnick_astrometry -t n1_2x2_slow_full.ecsv

  ``pipnick_reduce`` reprocesses only the images whose table values changed
  (see the reduction documentation), writing the corrected coordinates to
  their headers, and ``pipnick_astrometry`` then solves them with the usual
  matching.  This is the better choice if the reduced images will be used
  elsewhere, because it also corrects the coordinates in their headers.

Coordinate frames
-----------------

All solutions are in the ICRS, the frame of the Gaia catalog.

astrometry.net writes its solutions with ``EQUINOX = 2000.0`` and no
``RADESYS`` keyword, which by the FITS WCS standard means the FK5 (J2000)
frame, and its web API has no setting to request a different frame.
However, astrometry.net fits its solutions to reference catalogs (Tycho-2,
2MASS, and Gaia) whose positions are on the ICRS, so the solutions are
effectively in the ICRS; FK5 J2000 and the ICRS differ by at most a few tens
of milliarcseconds, far below the accuracy needed for Nickel data.
``pipnick`` therefore relabels astrometry.net solutions as ICRS
(``RADESYS = 'ICRS'``, without ``EQUINOX``) before applying them, which keeps
every solution in one frame.

Providing a solution
--------------------

If you know the pixel and celestial coordinates of a set of sources in an
image (at least 3), you can fit and apply a solution yourself.  Pixel
coordinates are 0-based (the numpy convention, as returned by ``photutils``).

.. code-block:: python

   from astropy.coordinates import SkyCoord
   from pipnick.astrometry.fit import fit_wcs
   from pipnick.pipelines.astrometry import apply_wcs

   coords = SkyCoord(ra, dec, unit='deg')
   wcs, stats = fit_wcs(x, y, coords)
   print(f'{stats["npts"]} sources used; RMS residual {stats["rms"]:.2f} arcsec')
   apply_wcs('rdx/M11_V_C006_2026-08-20.fits', wcs, stats=stats)

``fit_wcs`` fits a gnomonic (TAN) projection, optionally with SIP distortion
(``sip_degree``), and iteratively rejects outliers.  The applied solution is
marked ``PNASTMTH = 'points'``, is saved to the ``astrometry`` directory so
that it is restored if the image is rewritten, and is only replaced by
``pipnick_astrometry --resolve``.  Coordinates in the ICRS (e.g., from Gaia)
are preferred, but other frames are handled correctly.

Using the pipeline in Python
----------------------------

.. code-block:: python

   from pipnick.pipelines.astrometry import astrometry_all

   solved = astrometry_all('n1_2x2_slow_full.ecsv', max_submissions=1)

See :func:`~pipnick.pipelines.astrometry.astrometry_all` for the full set of
parameters.
