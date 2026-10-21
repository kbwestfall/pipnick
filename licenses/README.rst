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
  redistribute. Their licenses do not govern ``pipnick``'s own terms, but one
  of them constrains how ``pipnick`` may be packaged.


Adapted code
------------

PypeIt — ``PYPEIT_LICENSE.rst``
    BSD 3-Clause. https://github.com/pypeit/PypeIt

    Adapted into ``pipnick/pkg/logger.py``, ``pipnick/pkg/exceptions.py`` and
    ``pipnick/scripts/scriptbase.py``, and used as the template for this
    package's build and documentation configuration.

PypeIt permits redistribution with attribution, so the adapted code may be
shipped in the ``pipnick`` wheel; its license is included in the wheel's
metadata alongside ``LICENSE.rst``.


Dependencies
------------

.. warning::

    **loess is not redistributable, and no loess code may be copied into
    pipnick.**

``loess`` (Cappellari) is a dependency of ``pipnick`` — it supplies
``loess_2d``, used by ``pipnick.utils.graphs`` to smooth contour maps. Its
license is declared as ``Other/Proprietary License``, and reads in relevant
part:

    Copyright (C) 2010-2022, Michele Cappellari

    This software is provided as is without any warranty whatsoever.
    Permission to use, for non-commercial purposes is granted. Permission to
    modify for personal or internal use is granted, provided this copyright and
    disclaimer are included unchanged at the beginning of the file. All other
    rights are reserved. **In particular, redistribution of the code is not
    allowed.**

No copy of ``loess`` is included in this directory, because ``pipnick`` does
not redistribute it: ``pip`` installs it from PyPI as a declared dependency.
Two consequences follow:

1. **Never vendor, copy, or adapt loess source into pipnick.** Import it.
   Writing a clean-room implementation from the published method is
   permissible; a line-by-line translation is not.

2. **Expect packaging friction downstream.** conda-forge and distribution
   packagers cannot redistribute ``loess``, so anything that wants to ship
   ``pipnick`` as a system package will need the ``loess``-dependent path to
   be optional. ``loess`` is therefore imported only in
   ``pipnick/utils/graphs.py``.

All other dependencies (astropy, astroquery, ccdproc, IPython, matplotlib,
numpy, pandas, photutils, requests, scikit-image, and scipy) carry permissive
licenses — BSD, Apache-2.0, PSF, or equivalent — that impose no conditions on
``pipnick``.
