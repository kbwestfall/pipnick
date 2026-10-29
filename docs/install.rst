.. _install:

Installation
============

Requirements
------------

``pipnick`` requires Python 3.12, 3.13, or 3.14.  Its dependencies are
installed with it: ``numpy``, ``scipy``, ``astropy``, ``astroquery``,
``ccdproc``, ``photutils`` (3.1 or later), ``reproject``, ``scikit-image``,
``matplotlib``, and ``speclite``.

``pipnick`` is not yet released on PyPI; install it from GitHub.

Installing
----------

We recommend installing ``pipnick`` in its own virtual environment, kept
outside the directory of the source code:

.. code-block:: console

    python3 -m venv ~/.virtualenvs/pipnick
    source ~/.virtualenvs/pipnick/bin/activate

To install the latest version:

.. code-block:: console

    pip install "git+https://github.com/UCObservatories/pipnick.git"

To install from a clone of the repository, e.g., to follow its development or
to contribute:

.. code-block:: console

    git clone https://github.com/UCObservatories/pipnick.git
    cd pipnick
    pip install -e .

The optional dependencies are installed with ``pip install -e ".[test]"``
(running the tests), ``".[docs]"`` (building this documentation), or
``".[dev]"`` (both).

Checking the installation
-------------------------

The scripts are installed with the package.  Check that they run, and which
versions are installed:

.. code-block:: console

    pipnick_version --dependencies

Run the tests (with the ``test`` dependencies installed) with:

.. code-block:: console

    pytest --pyargs pipnick

The tests use synthetic data and never contact astrometry.net or the Gaia
archive.

.. _install-api-key:

The astrometry.net API key
--------------------------

The astrometric calibration (:doc:`astrometry`) submits some images to the
astrometry.net web service, which requires a free account.  Register at
https://nova.astrometry.net and copy your API key from your profile.
``pipnick`` reads the key from the ``ASTROMETRY_NET_API_KEY`` environment
variable, e.g., set in your shell startup file:

.. code-block:: console

    export ASTROMETRY_NET_API_KEY=<your key>

or, if that is not set, from the ``api_key`` entry in the ``[astrometry_net]``
section of your ``astroquery`` configuration file
(``~/.astropy/config/astroquery.cfg``):

.. code-block:: ini

    [astrometry_net]
    api_key = <your key>

The key is never a command-line argument, and is never logged.  The other
steps do not need it, and the astrometry only submits images that cannot be
solved otherwise.

Building the documentation
--------------------------

With the ``docs`` dependencies installed:

.. code-block:: console

    cd docs
    make html

The pages are written to ``docs/_build/html``.
