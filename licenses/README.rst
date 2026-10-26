Licenses
========

``pipnick`` is distributed under the BSD 3-Clause License; see ``LICENSE.rst``
in the repository root. This directory holds the licenses of third-party code
that ``pipnick`` either incorporates or depends on, together with the terms
each one imposes.

There are two distinct categories here, and the difference matters:

- **Adapted code** — source that has been copied into ``pipnick`` and modified.
  This is permitted only for permissively licensed projects, and every adapted
  module carries an attribution note in its module docstring naming the project
  and the file it came from.
- **Dependencies** — packages ``pipnick`` imports at run time but does not
  redistribute. Their licenses do not govern ``pipnick``'s own terms, but they
  can constrain how ``pipnick`` may be packaged.


Adapted code
------------

PypeIt — ``PYPEIT_LICENSE.rst``
    BSD 3-Clause. https://github.com/pypeit/PypeIt

    Adapted into ``pipnick/pkg/logger.py``, ``pipnick/pkg/exceptions.py`` and
    ``pipnick/scripts/scriptbase.py``, and used as the template for this
    package's build and documentation configuration.

astropy — ``ASTROPY_LICENSE.rst``
    BSD 3-Clause. https://github.com/astropy/astropy

    ``EllipticalMoffat2DModel`` in ``pipnick/photometry/moffat.py`` is adapted
    from ``Moffat2D`` in ``astropy/modeling/functional_models.py``.

MaNGA Data Analysis Pipeline (``mangadap``) — ``MANGADAP_LICENSE.rst``
    BSD 3-Clause. https://github.com/sdss/mangadap

    ``polygon_winding_number`` and ``point_inside_polygon`` in
    ``pipnick/utils/masking.py`` are adapted from ``mangadap/util/geometry.py``.

All three permit redistribution with attribution, so the adapted code may be
shipped in the ``pipnick`` wheel; their licenses are included in the wheel's
metadata alongside ``LICENSE.rst``.


Distributed data
----------------

PypeIt — ``PYPEIT_LICENSE.rst``
    ``pipnick/data/extinction/mthamextinct.dat`` (the Mt. Hamilton extinction
    curve) is copied unchanged from ``pypeit/data/extinction/mthamextinct.dat``.

Landolt standard stars
    ``pipnick/data/catalogs/landolt.ecsv`` is derived from the VizieR catalogs
    ``J/AJ/137/4186`` (Landolt 2009) and ``J/AJ/146/131`` (Landolt 2013),
    retrieved from the CDS.  The papers and VizieR are cited in the file's
    metadata and in ``pipnick.photometry.landolt``, with the acknowledgement
    requested by the CDS.


Dependencies
------------

All of the dependencies of ``pipnick`` (astropy, astroquery, ccdproc,
matplotlib, numpy, photutils, reproject, scikit-image, scipy, and speclite)
carry permissive licenses — BSD, Apache-2.0, PSF, or equivalent — that impose
no conditions on ``pipnick``.

.. warning::

    Do not add a dependency whose license restricts redistribution.  An
    earlier dependency, ``loess`` (Cappellari; used only by plotting code that
    has since been removed), forbids redistribution, which would have
    prevented ``pipnick`` from being packaged by conda-forge and Linux
    distributions.  No ``loess`` code may be copied into ``pipnick``.
