.. _astrometry:

Astrometric Calibration
=======================

After performing basic image processing, pipnick can construct an astrometric
solution for each science frame using ``pipnick_astrometry``.  These solutions
are attached to each science frame using a FITS World Coordinate System (WCS).
The code uses the `astrometry.net <https://nova.astrometry.net>`__ web service,
but submits as few images as possible: most images are solved locally, by
matching their sources to an image that is already solved or to the Gaia
catalog.

.. _astrometry-algorithm:

The algorithm
-------------

For each run of ``pipnick_astrometry``:

1. **Existing solutions are kept.**  Solutions that were lost because an image
   was rewritten (e.g., by rerunning ``pipnick_reduce``) are restored from the
   cache of earlier astrometry.net results or from solutions you provided
   (see :ref:`astrometry-user-solution`), without submitting anything.

2. **Images are grouped by pointing.**  Images whose telescope pointings are
   within ``--group_radius`` (default 1 arcmin) form a group.  A solution is
   found for one reference image per group, and the others are solved by
   matching their detected sources to those of the reference.  We refer to the
   latter as *bootstrapping* the solution.

3. Solutions for **groups without a solved image** are determined using an
   overlapping image with a solution, or by matching their sources to the Gaia
   DR3 catalog.  The match starts from a solution predicted using an existing
   solution from the same setup and the telescope pointing offset measured from
   existing solutions.

4. **Images that still cannot be solved are matched to Gaia with a wide
   search**, which allows for an incorrect telescope pointing in the
   metadata table (see :ref:`astrometry-incorrect-pointing`).

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

The results of the submissions are cached, keyed by what was submitted: the
positions of the sources, in order of brightness (or the image, with
``--upload_image``), and the pointing and pixel scale given to astrometry.net.
Any change to these changes the key, so that a cached result is no longer found
and the image could be submitted again: e.g., editing the pointing (``ra``,
``dec``) in the metadata table, or anything that changes the detected sources or
their ranking, such as ``--max_sources``, ``--detect_thresh``, or ``--fwhm``, but
also a version of ``photutils`` that measures the sources differently.  (``photutils`` 3.1 changed
the flux used to rank the sources, which is why pipnick requires it.)

Inputs
------

- **The ECSV metadata table** (``raw_table``), written by ``pipnick_metadata``;
  see :ref:`metadata`.  The pointing of each image (``ra``, ``dec``) is taken
  from it.

- **The reduced science images**, in the reduction directory.

- **An astrometry.net API key**, needed only if an image is submitted; see
  :ref:`install-api-key`.

- **Network access** to astrometry.net (for submissions) and to the Gaia archive
  (for queries, performed using ``astroquery``).  Saved results are reused, so a
  rerun usually needs neither.

- Optionally, **solutions you provide** for some images; see
  :ref:`astrometry-user-solution`.


Execution
---------

An example call looks like this:

.. code-block:: console

   pipnick_astrometry n1_2x2_slow_full.ecsv

The required argument is the per-setup metadata table.

The options are:

- ``-r``, ``--rdxdir``: The reduction directory.  If not provided, the
  directory saved in the metadata file is used.

- ``--dry_run``: Detect the sources, group the images, and report what would
  be submitted, without solving or submitting anything.  The source lists
  that would be submitted are written to ``astrometry/dry_run``.

- ``--max_submissions``: The maximum number of allowed submissions to
  astrometry.net.

- ``--images``: A limited set of images for which to determine astrometric
  solutions, given by their file names (with or without the extension) or
  combination group (``comb_id``).  The other images are not changed, but their
  solutions are used as references.

- ``--resolve``: Ignore the existing and cached solutions and solve the images
  (only those given by ``--images``, if provided) again.  This can submit
  images that were solved before.

- ``--no_qa``: Do not write the check figures to ``astrometry/qa``.  By
  default, a figure is redrawn for each image solved in the run, and drawn for
  any image without one (e.g., an unsolved image or every image in a first
  dry run).  It shows the image, the detected sources (circles, numbered by
  brightness), and, for an image solved in the run, the sources matched to the
  reference image or Gaia (crosses), with their residuals magnified 20 times
  (red).

- ``--max_sources``, ``--detect_thresh``, ``--fwhm``: The maximum number of
  (brightest) sources detected in each image (default 100), the detection
  threshold as a multiple of the background noise (default 5), and the FWHM
  of the sources in arcsec (default: typical for the camera).  Changing these
  changes the sources that would be submitted, so earlier astrometry.net
  results for an image are not reused.

- ``--upload_image``: Upload the images (with masked pixels set to 0) instead
  of source lists, so that astrometry.net detects the sources.

- ``--group_radius``, ``--min_matches``, ``--max_rms``: See
  :ref:`astrometry-algorithm`.

- ``--search``: The maximum error in arcsec for the expected position of an
  image when matching it to another image or to Gaia (default 60).

- ``--wide_search``, ``--wide_min_matches``: The maximum pointing error in
  arcmin allowed by the wide search (default 10; 0 disables it), and the
  minimum number of matched sources it requires (default 10).  See
  :ref:`astrometry-incorrect-pointing`.

A good first step for a new night is a dry run, followed by a run limited to
one submission:

.. code-block:: console

   pipnick_astrometry n1_2x2_slow_full.ecsv --dry_run
   pipnick_astrometry n1_2x2_slow_full.ecsv --max_submissions 1

Rerunning the script is safe: solved images are not changed, and nothing is
submitted unless an image cannot be solved otherwise.

Outputs
-------

The astrometric solution is written to the primary header of each reduced image,
replacing any existing WCS.  The following header keywords are also included and
record how the solution was found:

- ``PNASTMTH``: The method used: ``astrometry.net``, ``bootstrap`` (from another
  image), ``gaia``, or ``points`` (provided by you).

- ``PNASTWH``: A hash of the solution.

- ``PNASTVER``: The pipnick version that found the solution.

- ``PNASTNPT`` and ``PNASTRMS``: The number of matched sources and the RMS
  residual in arcsec (for bootstrapped, Gaia, and provided solutions).

- ``PNASTREF`` and ``PNASTRH``: The reference image (or ``Gaia DR3``) and the
  hash of its solution.  If the solution of a reference image changes, the
  images bootstrapped from it are solved again on the next run.

- ``PNASTIN`` and ``PNASTSUB``: The hash of the submitted input and the
  astrometry.net submission ID.

- ``PNASTPOE``: If the image center is farther than ``--search`` from the
  pointing in the metadata table, i.e., if the pointing is likely incorrect,
  the distance in arcmin between them.

The ``astrometry`` subdirectory of the reduction directory holds:

- ``anet_<hash>.fits``: Each solution returned by astrometry.net, saved
  exactly as returned (a FITS header), where ``<hash>`` identifies the
  submitted input.

- ``astrometry_net.ecsv``: An index of the astrometry.net submissions, with
  the input hash, image, status (``solved``, ``failed``, ``rejected``, or
  ``pending``), submission ID, solution file, and date of each.  Failed and
  rejected submissions are not resubmitted, and pending submissions are
  checked again on the next run.

- ``gaia_<n>.vot.gz``: Each Gaia query result, saved exactly as returned by
  the Gaia archive (a gzip-compressed VOTable).

- ``gaia_queries.ecsv``: An index of the Gaia queries, with the cone, the
  magnitude limit, the maximum number of sources, the date, and the query
  of each.  A saved query that covers an image is reused.

- ``<image>_points.fits``: The solutions you provided (see
  :ref:`astrometry-user-solution`).

- ``summary.ecsv``: The summary of the last run (see below).

- ``qa/<image>.png``: The quality assessment figures.

- ``dry_run/<image>_sources.ecsv``: The detected sources, written with
  ``--dry_run``.

The products returned by astrometry.net and Gaia are never modified; any
change (e.g., relabeling a solution as ICRS, see :ref:`astrometry-coo-frame`, or
propagating Gaia positions to the epoch of the observation) is applied when
they are used.  Do not delete this directory unless you want to submit the
images or query Gaia again.

The log and ``astrometry/summary.ecsv`` summarize the result for each image: its
status (``solved``, ``unsolved``, ``pending``, or ``skipped`` if not selected by
``--images``), method, reference, number of matched sources, RMS residual, and
pointing error.  For each unsolved image, the summary gives the reason it could
not be solved and a suggestion.  When the logging verbosity is set to ``-v 2``,
the log also gives the details of every attempt to match an image (the shift,
rotation, number of votes and matches, and RMS residual).

The log of each run is written to ``log/pipnick_astrometry_<time>.log`` in the
reduction directory; see :ref:`common`.

Troubleshooting
---------------

.. _astrometry-incorrect-pointing:

Incorrect telescope pointings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  area, at least ``--wide_min_matches`` sources (default 10) must match.
  Other images with the same (stale) pointing then use the measured error,
  so they are solved by the usual matching, even if they have too few sources
  for the wide search.

  Any image whose center is farther than ``--search`` (default 60 arcsec)
  from its pointing in the metadata table, however it was solved, is
  reported with a warning that gives the true image center, and the distance
  is recorded in the header (``PNASTPOE``).  These images are not used to
  measure the typical pointing offset used to predict the positions of the
  other images.

- **By correcting the metadata table.**  Edit the ``ra`` and ``dec`` of the
  affected frames in the metadata file (e.g., using the image centers
  reported by the wide search, or the coordinates in the observing log), and
  rerun the reduction and the astrometry:

  .. code-block:: console

     pipnick_reduce n1_2x2_slow_full.ecsv
     pipnick_astrometry n1_2x2_slow_full.ecsv

  ``pipnick_reduce`` reprocesses only the images whose table values changed
  (see :ref:`reduce`), writing the corrected coordinates to
  their headers, and ``pipnick_astrometry`` then solves them with the usual
  matching.  This is the better choice if the reduced images will be used
  elsewhere because it also corrects the coordinates in their headers.

Images that are not solved
~~~~~~~~~~~~~~~~~~~~~~~~~~

Finding astrometric solutions can be difficult, e.g., for sparse fields,
short exposures, or images with incorrect pointings.  If images are not
solved, or their solutions are poor:

- **Start with a dry run.**  ``--dry_run`` shows the pointing groups, the number
  of sources detected in each image, and which images would be submitted.  It
  also writes a figure of each image with its detected sources (unless the image
  already has one, or ``--no_qa`` is set), so you can check that they are stars
  and not, e.g., knots in a galaxy, cosmic rays, or artifacts along the image
  edges.  Images with fewer than about 10 sources (e.g., short exposures in the
  B band or sparse fields at high Galactic latitude) are the hardest to solve.

- **Read the summary.**  The reason each image was not solved, and a suggestion,
  are given at the end of the log and in ``astrometry/summary.ecsv``.

- **Check the pointing first.**  A warning that an image "is centered at ...
  arcmin from its pointing" usually means that the RA and Dec in the raw header
  are wrong; see :ref:`astrometry-incorrect-pointing`.  Chance matches rejected
  because of a large RMS residual are also a sign of a wrong pointing.

- **Allow for larger pointing errors.**  If images fail because too few of their
  sources match and the telescope pointing is poor, increase ``--search``
  (default 60 arcsec).  For larger, isolated errors, increase ``--wide_search``.

- **Get one good solution first.**  Every other method builds on a solved
  image in the same setup, so make sure the richest image is solved (e.g.,
  with ``--max_submissions 1``).  This is usually enough for the whole setup:
  images in other pointing groups are matched to Gaia using its solution.  If
  a pointing group still cannot be solved (e.g., in a sparse field), solving
  one image in that group (by a submission, or by providing a solution)
  allows the rest of the group to be bootstrapped from it.

- **Adjust the source detection.**  In sparse fields, a lower
  ``--detect_thresh`` finds more (fainter) sources; in crowded fields, a
  smaller ``--max_sources`` keeps only the brightest.  ``--fwhm`` should match
  the seeing.

- **Sparse fields.**  Lowering ``--min_matches`` makes false matches more
  likely.  Instead, provide a solution yourself from a handful of identified
  stars (see :ref:`astrometry-user-solution`), or let astrometry.net detect the
  sources with ``--upload_image``.

- **Rejected solutions.**  The log gives the reason, e.g., "RMS residual is
  0.41 arcsec".  In poor seeing or crowded fields, a modestly larger
  ``--max_rms`` (e.g., 0.5) may be acceptable for your science; check the
  image's figure in ``astrometry/qa`` first.

- **Iterate on the problem images.**  ``--images`` limits a run to the images
  you are working on.  Results are cached, and failed submissions are not
  repeated; ``--resolve`` ignores all of this and can submit again, so combine
  it with ``--images`` and ``--max_submissions``.

- **Inspect the result.**  Check ``PNASTMTH``, ``PNASTNPT``, ``PNASTRMS``, and
  ``PNASTPOE`` in the headers (or ``astrometry/summary.ecsv``), and the
  figures in ``astrometry/qa``.

.. _astrometry-coo-frame:

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

.. _astrometry-user-solution:

Providing a solution directly
-----------------------------

If you know the pixel and celestial coordinates of a set of sources in an
image (at least 3), you can fit and apply a solution yourself.  Pixel
coordinates are 0-based (the numpy convention, as returned by ``photutils``).

.. code-block:: python

   from astropy.coordinates import SkyCoord
   from pipnick.astrometry.fit import fit_wcs
   from pipnick.astrometry.header import apply_wcs

   coords = SkyCoord(ra, dec, unit='deg')
   wcs, stats = fit_wcs(x, y, coords)
   print(f'{stats["npts"]} sources used; RMS residual {stats["rms"]:.2f} arcsec')
   apply_wcs('rdx/M11_V_C006_2026-08-20.fits', wcs, stats=stats)

:func:`pipnick.astrometry.fit.fit_wcs` fits a gnomonic (TAN) projection,
optionally with SIP distortion (``sip_degree``), and iteratively rejects
outliers.  The applied solution is marked ``PNASTMTH = 'points'``, is saved to
the ``astrometry`` directory so that it is restored if the image is rewritten,
and is only replaced by ``pipnick_astrometry --resolve``.  Coordinates in the
ICRS (e.g., from Gaia) are preferred, but other frames are handled correctly.

Using the pipeline in Python
----------------------------

.. code-block:: python

   from pipnick.pipelines.astrometry import astrometry_all

   solved = astrometry_all('n1_2x2_slow_full.ecsv', max_submissions=1)

See :func:`~pipnick.pipelines.astrometry.astrometry_all` for the full set of
parameters.
